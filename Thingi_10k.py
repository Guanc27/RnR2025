import thingi10k

thingi10k.init() # Download the dataset and update cache

# Loop through all entries in the dataset
for entry in thingi10k.dataset():
    file_id = entry['file_id']
    author = entry['author']
    license = entry['licence']
    vertices, facets = thingi10k.load_file(entry['file_path'])
    # Do something with the vertices and facets

help(thingi10k) # for more information