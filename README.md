# zwaveformgen

This is more or less a Zig port of [cappelnord/waveformgen](https://github.com/cappelnord/waveformgen) that takes WAV audio data and generates:

* **PNG waveform visualizations** with customizable styling and rendering options
* **JSON amplitude data** for use in web applications or custom visualizations

## Features

- **Flexible waveform rendering**: Choose between min/max peaks or RMS (Root Mean Square) visualization
- **Customizable spacing**: Adjust horizontal spacing between waveform samples
- **Scalable output**: Configurable image height and amplitude scaling
- **Normalized RMS**: Option to normalize RMS values relative to the maximum amplitude for optimal dynamic range

## Installation

### Prerequisites

- [Zig](https://ziglang.org/) 0.14.1 or later

### Building from Source

```bash
zig build
```

The executable will be built to `zig-out/bin/zwaveformgen` (or `zwaveformgen.exe` on Windows).

## Usage

### Basic Usage

```bash
# Generate PNG waveform
./zwaveformgen --input audio.wav --output-type png --output-png waveform.png

# Generate JSON amplitude data
./zwaveformgen --input audio.wav --output-type json --output-json amplitudes.json
```

### Command Line Options

| Option | Description | Required | Default |
|--------|-------------|----------|---------|
| `--input <filename>` | Input WAV file | ✅ Yes | - |
| `--output-type <type>` | Output type: `png` or `json` | ✅ Yes | - |
| `--output-png <filename>` | Output PNG filename | ✅ Yes (if type=png) | - |
| `--output-json <filename>` | Output JSON filename | ✅ Yes (if type=json) | - |
| `--height <value>` | Image height in pixels | No | 280 |
| `--channel-spacing <value>` | Spacing between waveform samples | No | 1 |
| `--use-rms` | Use RMS values instead of min/max | No | false |
| `--normalize-rms` | Normalize RMS relative to maximum | No | false |
| `--help`, `-h` | Show help message | No | - |

## Output Formats

### PNG Images

- **Fixed width**: 1800 pixels (represents 1800 amplitude samples)
- **Configurable height**: Set with `--height` parameter (default: 280px)
- **RGBA format**: Transparent waveform on light gray background
- **Waveform rendering**:
  - **Min/Max mode** (default): Shows peak amplitudes as vertical lines
  - **RMS mode**: Shows Root Mean Square values as filled areas, optionally normalized

### JSON Output

The JSON output contains amplitude data suitable for web visualization:

```json
{
  "values": [23, 45, 12, 67, 89, ...],
  "count": 1800,
  "range": [0, 280]
}
```

- `values`: Array of 1800 amplitude values (0-255 range)
- `count`: Number of amplitude samples (always 1800)
- `range`: [min, max] values for the amplitude scale

## Audio Processing

- **Supported format**: 16-bit PCM WAV files
- **Channel handling**: Multi-channel audio is automatically mixed down
- **Downsampling**: Audio is processed to exactly 1800 amplitude points regardless of input length
- **RMS calculation**: Uses sliding window approach for smooth amplitude curves

## Waveform Generation Methods

1. **Min/Max (Default)**
2. **RMS (Root Mean Square)**

## Dependencies

- **zstbi**: Image processing library for PNG output
