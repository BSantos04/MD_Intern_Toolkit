from PIL import Image
import sys

# Define input parameters and open the input images
img1 = Image.open(sys.argv[1])
img2 = Image.open(sys.argv[2])
outfile = sys.argv[3]
img1 = img1.convert("RGBA")
img2 = img2.convert("RGBA")

# Make sure both images have the same size
if img1.size[1] != img2.size[1]:
   img1 = img1.resize((img1.size[0], img2.size[0]))

# Overlay both images in order to comapre them
new_img = Image.blend(img1,img2, 0.5)

# Save the images overlayed in a PNG file
new_img.save(outfile, "PNG")

