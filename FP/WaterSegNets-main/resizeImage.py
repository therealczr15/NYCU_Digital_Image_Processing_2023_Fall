import os, sys
from PIL import Image
from tqdm import tqdm

size = (256, 256)

image_path = 'training_dataset/mask/'
image_outpath = 'training_dataset/mask_reshape/'
image_names = sorted(next(os.walk(image_path))[-1])

for id in tqdm(range(len(image_names)), desc="Images"):
	path = image_path + image_names[id]
	outpath = image_outpath + image_names[id]
	try:
		im = Image.open(path)
		im = im.resize(size)
		im.save(outpath)
		
	except IOError:
		print("Error occured")