import os
from PIL import Image
import sys

fluid = 0b01000000
solid = 0b00100000
empty = 0b00010000

# Function to map value to color
def valToCol(value):
    if(value):
        return (255,255,255)
    else:
        return (0,0,0)

# Function to process each line and generate an image
def process_line(line, output_dir, idx):
    line = line.strip()  # Remove whitespace or newline characters
    if not line or ';' not in line:
        return  # Ignore empty or malformed lines

    # Parse the line
    try:
        size, contents = line.split(':')
        width, height = map(int, size.split('x'))
        values = list(map(int, contents.strip(';').split(',')))
        
        if len(values) != width * height:
            raise ValueError("Grid size does not match the number of values")
    except Exception as e:
        print(f"Error parsing line: {line}\n{e}")
        return

    # Create an image
    img = Image.new('RGB', (width, height))
    pixels = img.load()

    for y in range(height):
        for x in range(width):
            value = values[y * width + x]
            pixels[x, y] = valToCol(value)
            #print(f"Value: {value}, color: {valToCol(value)}")

    # Save the image
    output_file = os.path.join(output_dir, f"grid_{width}x{height}_{idx}.png")
    img.save(output_file)
    print(f"Saved image: {output_file}")

# Main function to process the file
def main(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    i = 0
    with open(input_file, 'r') as file:
        for line in file:
            process_line(line, output_dir, i)
            i+=1

# Specify input file and output directory
input_file = sys.argv[1]
output_dir = './'

main(input_file, output_dir)
