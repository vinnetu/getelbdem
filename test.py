import getelbdem

fn = "tartu.tif"
x1 = 657500
y1 = 6473500
x2 = 661000
y2 = 6475700


# Download the file to wanted location defined by L-EST97 coordinate system
# Resolution can be either 1, 5, or 10
getelbdem.dlarea(x1, y1, x2, y2, fn, resolution = 5)