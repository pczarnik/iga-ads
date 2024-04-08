from PIL import Image
img_rgb = Image.open('forest_map.png').convert('L')
img_rgb.save('forest_map_grey.png')
