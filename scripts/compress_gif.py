#!/usr/bin/env python3
"""
GIF Compression Utility
Reduces GIF file size using various optimization techniques.

Usage:
    python compress_gif.py input.gif [output.gif] [--quality QUALITY] [--colors COLORS] [--scale SCALE]

Options:
    --quality QUALITY   Compression quality (1-3): 1=best quality, 3=smallest size (default: 2)
    --colors COLORS     Number of colors to use (8-256, default: 128)
    --scale SCALE       Scale factor (0.1-1.0, default: 1.0 for no scaling)
    --fps FPS          Target frames per second (default: auto-detect)
    --skip SKIP        Skip every N frames (default: 1, keep all frames)
"""

import sys
import os
import argparse
from pathlib import Path
from PIL import Image, ImageSequence

def get_file_size_mb(filepath):
    """Get file size in MB."""
    size_bytes = os.path.getsize(filepath)
    return size_bytes / (1024 * 1024)

def compress_gif(input_path, output_path=None, quality=2, colors=128, scale=1.0, skip_frames=1, target_fps=None):
    """
    Compress a GIF file.
    
    Args:
        input_path: Path to input GIF
        output_path: Path to output GIF (default: input_compressed.gif)
        quality: Compression quality (1=best, 3=smallest)
        colors: Number of colors (8-256)
        scale: Scale factor (0.1-1.0)
        skip_frames: Skip every N frames (1=keep all)
        target_fps: Target frames per second (None=keep original)
    """
    input_path = Path(input_path)
    
    if not input_path.exists():
        print(f"❌ Error: File not found: {input_path}")
        return False
    
    if output_path is None:
        output_path = input_path.parent / f"{input_path.stem}_compressed.gif"
    else:
        output_path = Path(output_path)
    
    print("=" * 70)
    print("                    GIF Compression Utility")
    print("=" * 70)
    print(f"Input:       {input_path}")
    print(f"Output:      {output_path}")
    print(f"Quality:     {quality} (1=best, 3=smallest)")
    print(f"Colors:      {colors}")
    print(f"Scale:       {scale:.2f}")
    print(f"Skip frames: {skip_frames}")
    print()
    
    # Load the GIF
    print("Loading GIF...")
    img = Image.open(input_path)
    
    # Get original properties
    original_size = get_file_size_mb(input_path)
    try:
        original_duration = img.info.get('duration', 100)
    except:
        original_duration = 100
    
    original_frames = 0
    for _ in ImageSequence.Iterator(img):
        original_frames += 1
    
    print(f"Original size:   {original_size:.2f} MB")
    print(f"Original frames: {original_frames}")
    print(f"Original size:   {img.size[0]}x{img.size[1]}")
    print()
    
    # Process frames
    print("Processing frames...")
    img.seek(0)
    frames = []
    durations = []
    frame_count = 0
    
    for i, frame in enumerate(ImageSequence.Iterator(img)):
        # Skip frames if requested
        if i % skip_frames != 0:
            continue
        
        # Convert to RGB mode if needed
        if frame.mode != 'RGB':
            frame = frame.convert('RGB')
        
        # Scale frame if requested
        if scale != 1.0:
            new_size = (int(frame.size[0] * scale), int(frame.size[1] * scale))
            frame = frame.resize(new_size, Image.Resampling.LANCZOS)
        
        # Reduce colors using adaptive palette
        frame = frame.convert('P', palette=Image.Palette.ADAPTIVE, colors=colors)
        
        frames.append(frame)
        
        # Calculate duration
        if target_fps:
            duration = int(1000 / target_fps)
        else:
            duration = frame.info.get('duration', original_duration)
        durations.append(duration)
        
        frame_count += 1
        if frame_count % 10 == 0:
            print(f"  Processed {frame_count} frames...", end='\r')
    
    print(f"  Processed {frame_count} frames... Done!")
    print()
    
    # Save compressed GIF
    print("Saving compressed GIF...")
    
    # Set optimization level based on quality
    optimize = True
    if quality == 1:
        # Best quality
        save_kwargs = {'optimize': optimize, 'quality': 95}
    elif quality == 2:
        # Balanced
        save_kwargs = {'optimize': optimize, 'quality': 85}
    else:
        # Smallest size
        save_kwargs = {'optimize': optimize, 'quality': 75}
    
    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        **save_kwargs
    )
    
    # Report results
    compressed_size = get_file_size_mb(output_path)
    reduction = ((original_size - compressed_size) / original_size) * 100
    
    print()
    print("=" * 70)
    print("                         Results")
    print("=" * 70)
    print(f"Original:    {original_size:.2f} MB ({original_frames} frames)")
    print(f"Compressed:  {compressed_size:.2f} MB ({frame_count} frames)")
    print(f"Reduction:   {reduction:.1f}%")
    print(f"Saved:       {original_size - compressed_size:.2f} MB")
    print("=" * 70)
    print(f"✓ Output saved to: {output_path}")
    print()
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Compress GIF files to reduce size',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic compression (default settings)
  python compress_gif.py animation.gif
  
  # High quality compression with 256 colors
  python compress_gif.py animation.gif output.gif --quality 1 --colors 256
  
  # Maximum compression
  python compress_gif.py animation.gif output.gif --quality 3 --colors 64 --scale 0.75
  
  # Skip every other frame and reduce to 15 fps
  python compress_gif.py animation.gif output.gif --skip 2 --fps 15
        """
    )
    
    parser.add_argument('input', help='Input GIF file')
    parser.add_argument('output', nargs='?', help='Output GIF file (default: input_compressed.gif)')
    parser.add_argument('--quality', type=int, default=2, choices=[1, 2, 3],
                        help='Compression quality: 1=best quality, 2=balanced, 3=smallest size (default: 2)')
    parser.add_argument('--colors', type=int, default=128,
                        help='Number of colors to use (8-256, default: 128)')
    parser.add_argument('--scale', type=float, default=1.0,
                        help='Scale factor (0.1-1.0, default: 1.0)')
    parser.add_argument('--fps', type=int, default=None,
                        help='Target frames per second (default: keep original)')
    parser.add_argument('--skip', type=int, default=1,
                        help='Skip every N frames (default: 1, keep all frames)')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.colors < 8 or args.colors > 256:
        print("❌ Error: Colors must be between 8 and 256")
        return 1
    
    if args.scale < 0.1 or args.scale > 1.0:
        print("❌ Error: Scale must be between 0.1 and 1.0")
        return 1
    
    if args.skip < 1:
        print("❌ Error: Skip must be at least 1")
        return 1
    
    # Compress the GIF
    success = compress_gif(
        args.input,
        args.output,
        quality=args.quality,
        colors=args.colors,
        scale=args.scale,
        skip_frames=args.skip,
        target_fps=args.fps
    )
    
    return 0 if success else 1

if __name__ == '__main__':
    sys.exit(main())
