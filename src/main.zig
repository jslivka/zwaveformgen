const std = @import("std");
const zstbi = @import("zstbi");

const OutputType = enum {
    png,
    json,
};

const ChunkHeader = struct {
    id: [4]u8,
    size: u32,

    pub fn read(reader: anytype) !ChunkHeader {
        var header: ChunkHeader = undefined;
        try reader.readNoEof(&header.id);
        header.size = try reader.readInt(u32, .little);
        return header;
    }
};

const WavHeader = struct {
    riff_header: ChunkHeader,
    wave_format: [4]u8,

    pub fn read(reader: anytype) !WavHeader {
        var header: WavHeader = undefined;
        header.riff_header = try ChunkHeader.read(reader);
        try reader.readNoEof(&header.wave_format);
        return header;
    }

    pub fn validate(self: WavHeader) !void {
        if (!std.mem.eql(u8, &self.riff_header.id, &[_]u8{ 'R', 'I', 'F', 'F' }) or
            !std.mem.eql(u8, &self.wave_format, &[_]u8{ 'W', 'A', 'V', 'E' }))
        {
            return error.InvalidWavFormat;
        }
    }
};

const FormatChunk = struct {
    header: ChunkHeader,
    audio_format: u16,
    num_channels: u16,
    sample_rate: u32,
    byte_rate: u32,
    block_align: u16,
    bits_per_sample: u16,

    pub fn read(reader: anytype) !FormatChunk {
        var chunk: FormatChunk = undefined;
        chunk.header = try ChunkHeader.read(reader);
        if (!std.mem.eql(u8, &chunk.header.id, &[_]u8{ 'f', 'm', 't', ' ' })) {
            return error.UnexpectedChunk;
        }
        chunk.audio_format = try reader.readInt(u16, .little);
        chunk.num_channels = try reader.readInt(u16, .little);
        chunk.sample_rate = try reader.readInt(u32, .little);
        chunk.byte_rate = try reader.readInt(u32, .little);
        chunk.block_align = try reader.readInt(u16, .little);
        chunk.bits_per_sample = try reader.readInt(u16, .little);

        return chunk;
    }
};

const DataChunk = struct {
    header: ChunkHeader,

    pub fn readSamples(self: DataChunk, reader: anytype, fmt_chunk: FormatChunk, allocator: std.mem.Allocator) ![]u8 {
        if (fmt_chunk.num_channels == 0 or fmt_chunk.bits_per_sample % 8 != 0) {
            return error.InvalidAudioFormat;
        }

        // Allocate buffer for all samples
        const buffer = try allocator.alloc(u8, self.header.size);
        errdefer allocator.free(buffer);

        const bytes_read = try reader.readAll(buffer);
        if (bytes_read != self.header.size) {
            return error.UnexpectedEOF;
        }

        return buffer;
    }

    pub fn read(reader: anytype) !DataChunk {
        var chunk: DataChunk = undefined;
        chunk.header = try ChunkHeader.read(reader);
        if (!std.mem.eql(u8, &chunk.header.id, &[_]u8{ 'd', 'a', 't', 'a' })) {
            return error.UnexpectedChunk;
        }
        return chunk;
    }
};

fn jsonStrFromSamplesArray(samples_array: []u8, image_height: usize, allocator: std.mem.Allocator) ![]u8 {
    var string = std.ArrayList(u8).init(allocator);
    const writer = string.writer();

    // Create compact JSON object
    try writer.writeAll("{\"values\":[");

    for (samples_array, 0..) |sample, i| {
        if (i > 0) {
            try writer.writeAll(",");
        }
        try writer.print("{}", .{sample});
    }

    try writer.writeAll("],\"count\":");
    try writer.print("{}", .{samples_array.len});
    try writer.writeAll(",\"range\":[0,");
    // Report the actual effective range (clamped to u8 max)
    const effective_range = @min(image_height, 255);
    try writer.print("{}", .{effective_range});
    try writer.writeAll("]}");
    return string.toOwnedSlice();
}

// Function to generate an array of amplitude values from WAV samples
// This function takes raw WAV sample data and downsamples it to exactly 1800 amplitude values.
// Each value represents the root mean square (RMS) amplitude normalized relative to the maximum RMS found in the waveform.
// Values are scaled to the specified image height but clamped to u8 range (0-255) for storage.
//
// Parameters:
//   samples: Raw byte data from the WAV file's data chunk
//   fmt_chunk: Format information from the WAV file (sample rate, bit depth, etc.)
//   image_height: Target height for scaling (values will be clamped to 0-255 range)
//   allocator: Memory allocator for the output array
//
// Returns: An allocated array of 1800 u8 values representing normalized amplitudes in range [0, min(image_height, 255)]
//          Caller is responsible for freeing this memory with allocator.free()
pub fn generateArrayFromSamples(samples: []const u8, fmt_chunk: FormatChunk, image_height: usize, allocator: std.mem.Allocator) ![]u8 {
    const target_count = 1800;
    const max_amplitude_value = image_height;

    // Allocate output array
    var amplitudes = try allocator.alloc(u8, target_count);
    errdefer allocator.free(amplitudes);

    // Allocate temporary array to store raw RMS values
    var rms_values = try allocator.alloc(f64, target_count);
    defer allocator.free(rms_values);

    // Calculate samples per pixel
    const bytes_per_sample = fmt_chunk.bits_per_sample / 8;
    const samples_per_channel = samples.len / (fmt_chunk.num_channels * bytes_per_sample);
    const samples_per_output = samples_per_channel / target_count;

    if (samples_per_output == 0) return error.InsufficientSamples;

    // First pass: Calculate all RMS values and find maximum
    var max_rms: f64 = 0.0;
    var i: usize = 0;
    while (i < target_count) : (i += 1) {
        var sum_squares: f64 = 0.0;

        // Sum the squares of sample values in this output element's sample range
        var s: usize = 0;
        while (s < samples_per_output) : (s += 1) {
            const sample_idx = (i * samples_per_output + s) * fmt_chunk.num_channels * bytes_per_sample;
            if (sample_idx + 1 >= samples.len) break;

            // Convert bytes to 16-bit sample (assuming PCM)
            const sample_value = @as(i16, @bitCast([2]u8{
                samples[sample_idx],
                samples[sample_idx + 1],
            }));

            const sample_f64: f64 = @floatFromInt(sample_value);
            sum_squares += sample_f64 * sample_f64;
        }

        // Calculate RMS as the square root of the average of squared values
        const rms = @sqrt(sum_squares / @as(f64, @floatFromInt(samples_per_output)));
        rms_values[i] = rms;

        // Track maximum RMS value
        if (rms > max_rms) {
            max_rms = rms;
        }
    }

    // Second pass: Normalize relative to maximum RMS found
    i = 0;
    while (i < target_count) : (i += 1) {
        if (max_rms > 0.0) {
            // Normalize to 0-1 range relative to the maximum RMS found
            const normalized = rms_values[i] / max_rms;

            // Scale to image_height but clamp to u8 range for storage
            const max_amplitude_f64 = @as(f64, @floatFromInt(max_amplitude_value));
            const scaled_value = normalized * max_amplitude_f64;

            // Clamp to u8 range (0-255) to prevent overflow
            const clamped_value = @min(@max(scaled_value, 0.0), 255.0);
            amplitudes[i] = @as(u8, @intFromFloat(clamped_value));
        } else {
            // Handle case where all samples are silent
            amplitudes[i] = 0;
        }
    }

    return amplitudes;
}

const WaveformImage = struct {
    width: usize,
    height: usize,
    channel_spacing: usize = 1, // Spacing between samples (1 = no spacing, 2 = every other pixel, etc.)
    // Higher values create more spaced out waveforms with gaps between sample lines
    use_rms: bool = false, // Whether to use RMS values for waveform height instead of min/max
    normalize_rms: bool = false, // Whether to normalize RMS values relative to maximum RMS found
    waveform_color: [4]u8 = .{ 0, 0, 0, 0 }, // Transparent
    background_color: [4]u8 = .{ 239, 239, 239, 255 }, // Offwhite

    pub fn generateFromSamples(self: WaveformImage, samples: []const u8, fmt_chunk: FormatChunk, output_filename: []const u8, allocator: std.mem.Allocator) !void {
        // Create image buffer (RGBA format)
        const image_size = self.width * self.height * 4;
        var image = try allocator.alloc(u8, image_size);
        defer allocator.free(image);

        // Fill with background color
        var i: usize = 0;
        while (i < image_size) : (i += 4) {
            image[i] = self.background_color[0];
            image[i + 1] = self.background_color[1];
            image[i + 2] = self.background_color[2];
            image[i + 3] = self.background_color[3];
        }

        // Calculate samples per pixel, accounting for channel spacing
        const bytes_per_sample = fmt_chunk.bits_per_sample / 8;
        const samples_per_channel = samples.len / (fmt_chunk.num_channels * bytes_per_sample);

        // Calculate effective width (number of sample groups to display)
        const effective_samples = self.width / self.channel_spacing;
        const samples_per_pixel = samples_per_channel / effective_samples;
        if (samples_per_pixel == 0) return error.ImageTooWide;

        // Process samples and draw waveform with spacing
        var sample_idx: usize = 0;
        while (sample_idx < effective_samples) : (sample_idx += 1) {
            const x = sample_idx * self.channel_spacing;
            if (x >= self.width) break;

            if (self.use_rms) {
                if (self.normalize_rms) {
                    // Two-pass approach for normalized RMS visualization
                    // First pass: Calculate all RMS values and find maximum
                    var rms_values = try allocator.alloc(f64, effective_samples);
                    defer allocator.free(rms_values);
                    var max_rms: f64 = 0.0;

                    var pass1_idx: usize = 0;
                    while (pass1_idx < effective_samples) : (pass1_idx += 1) {
                        var sum_squares: f64 = 0.0;
                        var sample_count: usize = 0;

                        var s: usize = 0;
                        while (s < samples_per_pixel) : (s += 1) {
                            const sample_byte_idx = (pass1_idx * samples_per_pixel + s) * fmt_chunk.num_channels * bytes_per_sample;
                            if (sample_byte_idx + 1 >= samples.len) break;

                            const sample_value = @as(i16, @bitCast([2]u8{
                                samples[sample_byte_idx],
                                samples[sample_byte_idx + 1],
                            }));

                            const sample_f64: f64 = @floatFromInt(sample_value);
                            sum_squares += sample_f64 * sample_f64;
                            sample_count += 1;
                        }

                        if (sample_count > 0) {
                            const rms = @sqrt(sum_squares / @as(f64, @floatFromInt(sample_count)));
                            rms_values[pass1_idx] = rms;
                            if (rms > max_rms) {
                                max_rms = rms;
                            }
                        } else {
                            rms_values[pass1_idx] = 0.0;
                        }
                    }

                    // Second pass: Draw normalized RMS waveform
                    sample_idx = 0;
                    while (sample_idx < effective_samples) : (sample_idx += 1) {
                        const x_pos = sample_idx * self.channel_spacing;
                        if (x_pos >= self.width) break;

                        if (max_rms > 0.0) {
                            const normalized = rms_values[sample_idx] / max_rms;
                            const rms_height = @as(usize, @intFromFloat(normalized * @as(f64, @floatFromInt(self.height))));

                            // Draw normalized RMS-based waveform (centered vertically)
                            const center_y = self.height / 2;
                            const half_height = @min(rms_height / 2, center_y);

                            var y: usize = center_y - half_height;
                            const end_y = center_y + half_height;
                            while (y <= end_y and y < self.height) : (y += 1) {
                                const pixel_idx = (y * self.width + x_pos) * 4;
                                image[pixel_idx] = self.waveform_color[0];
                                image[pixel_idx + 1] = self.waveform_color[1];
                                image[pixel_idx + 2] = self.waveform_color[2];
                                image[pixel_idx + 3] = self.waveform_color[3];
                            }
                        }
                    }
                } else {
                    // Original single-pass RMS approach
                    while (sample_idx < effective_samples) : (sample_idx += 1) {
                        const x_pos = sample_idx * self.channel_spacing;
                        if (x_pos >= self.width) break;

                        // Calculate RMS for this sample group
                        var sum_squares: f64 = 0.0;
                        var sample_count: usize = 0;

                        var s: usize = 0;
                        while (s < samples_per_pixel) : (s += 1) {
                            const sample_byte_idx = (sample_idx * samples_per_pixel + s) * fmt_chunk.num_channels * bytes_per_sample;
                            if (sample_byte_idx + 1 >= samples.len) break;

                            // Convert bytes to 16-bit sample (assuming PCM)
                            const sample_value = @as(i16, @bitCast([2]u8{
                                samples[sample_byte_idx],
                                samples[sample_byte_idx + 1],
                            }));

                            const sample_f64: f64 = @floatFromInt(sample_value);
                            sum_squares += sample_f64 * sample_f64;
                            sample_count += 1;
                        }

                        if (sample_count > 0) {
                            // Calculate RMS and scale to image height
                            const rms = @sqrt(sum_squares / @as(f64, @floatFromInt(sample_count)));
                            const normalized = rms / 32767.0; // Normalize to 0-1 range
                            const rms_height = @as(usize, @intFromFloat(normalized * @as(f64, @floatFromInt(self.height))));

                            // Draw RMS-based waveform (centered vertically)
                            const center_y = self.height / 2;
                            const half_height = @min(rms_height / 2, center_y);

                            var y: usize = center_y - half_height;
                            const end_y = center_y + half_height;
                            while (y <= end_y and y < self.height) : (y += 1) {
                                const pixel_idx = (y * self.width + x_pos) * 4;
                                image[pixel_idx] = self.waveform_color[0];
                                image[pixel_idx + 1] = self.waveform_color[1];
                                image[pixel_idx + 2] = self.waveform_color[2];
                                image[pixel_idx + 3] = self.waveform_color[3];
                            }
                        }
                    }
                }
            } else {
                // Original min/max approach
                var min: i16 = 32767;
                var max: i16 = -32768;

                // Find min/max values for this pixel's sample range
                var s: usize = 0;
                while (s < samples_per_pixel) : (s += 1) {
                    const sample_byte_idx = (sample_idx * samples_per_pixel + s) * fmt_chunk.num_channels * bytes_per_sample;
                    if (sample_byte_idx + 1 >= samples.len) break;

                    // Convert bytes to 16-bit sample (assuming PCM)
                    const sample_value = @as(i16, @bitCast([2]u8{
                        samples[sample_byte_idx],
                        samples[sample_byte_idx + 1],
                    }));

                    min = @min(min, sample_value);
                    max = @max(max, sample_value);
                }

                // Scale to fit height
                const sample_range = 65536; // Total range of 16-bit audio (-32768 to 32767)
                const height_i64 = @as(i64, @intCast(self.height));
                const scaled_min = @as(usize, @intCast(@divTrunc((@as(i64, min) + 32768) * height_i64, sample_range)));
                const scaled_max = @as(usize, @intCast(@divTrunc((@as(i64, max) + 32768) * height_i64, sample_range)));

                // Draw vertical line
                var y: usize = scaled_min;
                while (y <= scaled_max) : (y += 1) {
                    const pixel_idx = (y * self.width + x) * 4;
                    image[pixel_idx] = self.waveform_color[0];
                    image[pixel_idx + 1] = self.waveform_color[1];
                    image[pixel_idx + 2] = self.waveform_color[2];
                    image[pixel_idx + 3] = self.waveform_color[3];
                }
            }
        }

        // Create a zstbi.Image from raw buffer
        var img = zstbi.Image{
            .data = image,
            .width = @intCast(self.width),
            .height = @intCast(self.height),
            .num_components = 4,
            .bytes_per_component = 1,
            .bytes_per_row = @intCast(self.width * 4),
            .is_hdr = false,
        };
        // Save the image
        const output_filename_z = try allocator.dupeZ(u8, output_filename);
        defer allocator.free(output_filename_z);
        try img.writeToFile(output_filename_z, .png);
        std.debug.print("\nWaveform image saved as '{s}'\n", .{output_filename});
    }
};

pub fn main() !void {
    // Initialize allocator and zstbi
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse command line arguments
    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    var channel_spacing: usize = 1; // Default: no spacing
    var use_rms: bool = false; // Default: use min/max approach
    var normalize_rms: bool = false; // Default: don't normalize RMS
    var image_height: usize = 280; // Default height
    var output_png: []const u8 = ""; // Required: must specify output PNG file
    var output_json: []const u8 = ""; // Required: must specify output JSON file
    var output_type: ?OutputType = null; // Required: must specify output type
    var input_file: []const u8 = ""; // Required: must specify input file

    // Check for help flag
    for (args) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            std.debug.print("Usage: zwaveformgen [options]\n", .{});
            std.debug.print("Options:\n", .{});
            std.debug.print("  --input <filename>         Set input WAV file (required)\n", .{});
            std.debug.print("  --channel-spacing <value>  Set spacing between waveform samples (default: 1)\n", .{});
            std.debug.print("                             Higher values create more spaced out waveforms\n", .{});
            std.debug.print("  --height <value>           Set image height in pixels (default: 280)\n", .{});
            std.debug.print("                             RMS values will be scaled to this height\n", .{});
            std.debug.print("  --output-type <type>       Set output type: 'png' or 'json' (required)\n", .{});
            std.debug.print("  --output-png <filename>    Set output PNG filename (required if output is png)\n", .{});
            std.debug.print("  --output-json <filename>   Set output JSON filename (required if output is json)\n", .{});
            std.debug.print("  --use-rms                  Use RMS values for waveform height (default: min/max)\n", .{});
            std.debug.print("                             RMS values are scaled relative to image height\n", .{});
            std.debug.print("  --normalize-rms            Normalize RMS values relative to maximum RMS found\n", .{});
            std.debug.print("                             Ensures loudest part uses full dynamic range\n", .{});
            std.debug.print("  --help, -h                 Show this help message\n", .{});
            return;
        }
    }

    // Parse arguments
    for (args, 0..) |arg, i| {
        if (std.mem.eql(u8, arg, "--input") and i + 1 < args.len) {
            input_file = args[i + 1];
        } else if (std.mem.eql(u8, arg, "--use-rms")) {
            use_rms = true;
        } else if (std.mem.eql(u8, arg, "--normalize-rms")) {
            normalize_rms = true;
        } else if (std.mem.eql(u8, arg, "--height") and i + 1 < args.len) {
            image_height = std.fmt.parseInt(usize, args[i + 1], 10) catch blk: {
                std.debug.print("Invalid height value. Using default value of 280.\n", .{});
                break :blk 280;
            };
            if (image_height == 0) {
                std.debug.print("Height must be greater than 0. Using default value of 280.\n", .{});
                image_height = 280;
            }
        } else if (std.mem.eql(u8, arg, "--output-png") and i + 1 < args.len) {
            output_png = args[i + 1];
        } else if (std.mem.eql(u8, arg, "--output-json") and i + 1 < args.len) {
            output_json = args[i + 1];
        } else if (std.mem.eql(u8, arg, "--output-type") and i + 1 < args.len) {
            const type_str = args[i + 1];
            if (std.mem.eql(u8, type_str, "png")) {
                output_type = .png;
            } else if (std.mem.eql(u8, type_str, "json")) {
                output_type = .json;
            } else {
                std.debug.print("Invalid output type '{s}'. Valid options are: 'png', 'json'.\n", .{type_str});
                return;
            }
        } else if (std.mem.eql(u8, arg, "--channel-spacing") and i + 1 < args.len) {
            channel_spacing = std.fmt.parseInt(usize, args[i + 1], 10) catch blk: {
                std.debug.print("Invalid channel spacing value. Using default value of 1.\n", .{});
                break :blk 1;
            };
            if (channel_spacing == 0) {
                std.debug.print("Channel spacing must be greater than 0. Using default value of 1.\n", .{});
                channel_spacing = 1;
            }
            break;
        }
    }

    // Validate required arguments
    if (output_type == null) {
        std.debug.print("Error: --output-type is required. Use 'png' or 'json'.\n", .{});
        std.debug.print("Run with --help for usage information.\n", .{});
        return;
    }

    if (input_file.len == 0) {
        std.debug.print("Error: --input is required. Specify the input WAV file.\n", .{});
        std.debug.print("Run with --help for usage information.\n", .{});
        return;
    }

    // Validate that the corresponding output file parameter is provided
    if (output_type.? == .png and output_png.len == 0) {
        std.debug.print("Error: --output-png is required when --output-type is 'png'.\n", .{});
        std.debug.print("Run with --help for usage information.\n", .{});
        return;
    }

    if (output_type.? == .json and output_json.len == 0) {
        std.debug.print("Error: --output-json is required when --output-type is 'json'.\n", .{});
        std.debug.print("Run with --help for usage information.\n", .{});
        return;
    }

    if (channel_spacing > 1) {
        std.debug.print("Using channel spacing: {}\n", .{channel_spacing});
    }
    if (use_rms) {
        std.debug.print("Using RMS values for waveform height scaling\n", .{});
    }
    if (normalize_rms) {
        std.debug.print("Using normalized RMS values relative to maximum amplitude\n", .{});
    }
    if (image_height != 280) {
        std.debug.print("Using image height: {}\n", .{image_height});
    }

    zstbi.init(allocator);
    defer zstbi.deinit();

    // Open the audio file
    var audio_file = try std.fs.cwd().openFile(input_file, .{});
    defer audio_file.close();

    const reader = audio_file.reader();

    // Read and validate WAV header
    const wav_header = try WavHeader.read(reader);
    try wav_header.validate();

    // Read format chunk
    const format_chunk = try FormatChunk.read(reader);

    // Read chunks until we find the data chunk and generate waveform
    var audio_data: []u8 = undefined;
    defer allocator.free(audio_data);

    while (true) {
        const chunk_header = try ChunkHeader.read(reader);
        if (std.mem.eql(u8, &chunk_header.id, &[_]u8{ 'd', 'a', 't', 'a' })) {
            const data_chunk = DataChunk{ .header = chunk_header };
            audio_data = try data_chunk.readSamples(reader, format_chunk, allocator);
            break;
        }
        // Skip other chunks
        try reader.skipBytes(chunk_header.size, .{});
    }

    // Generate outputs based on selected type
    if (output_type.? == .png) {
        // Generate waveform visualization with configurable options
        var waveform = WaveformImage{
            .width = 1800,
            .height = image_height,
            .channel_spacing = channel_spacing,
            .use_rms = use_rms,
            .normalize_rms = normalize_rms,
        };
        try waveform.generateFromSamples(audio_data, format_chunk, output_png, allocator);
    } else if (output_type.? == .json) {
        // Generate amplitude array from samples (scaled to image height)
        const amplitudes = try generateArrayFromSamples(audio_data, format_chunk, image_height, allocator);
        defer allocator.free(amplitudes);

        const amplitudes_json = try jsonStrFromSamplesArray(amplitudes, image_height, allocator);
        defer allocator.free(amplitudes_json);

        // Write amplitudes to JSON file
        const output_json_z = try allocator.dupeZ(u8, output_json);
        defer allocator.free(output_json_z);
        const file = try std.fs.cwd().createFile(output_json_z, .{});
        defer file.close();

        try file.writeAll(amplitudes_json);
        std.debug.print("Amplitude values saved to '{s}'\n", .{output_json});
    }
}
