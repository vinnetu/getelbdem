from scipy.spatial import distance
import os.path
import subprocess
import fiona
import math
import copy
import rasterio
from rasterio.mask import mask
from shapely.geometry import box
from fiona.crs import from_epsg
import geopandas as gpd
import json
import requests
import os

def cropraster(fp: str, minx: float, miny: float, maxx: float, maxy: float, epsg_code: float = 3301) -> bool:
    # cuts a raster from big raster with filename generated from site ID

    data = rasterio.open(fp)
    # Create bounding box
    bbox = box(minx, miny, maxx, maxy)

    geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=from_epsg(epsg_code))

    coords = [json.loads(geo.to_json())['features'][0]['geometry']]
    out_img, out_transform = mask(dataset=data, shapes=coords, crop=True)
    out_meta = data.meta.copy()
    print(data)
    out_meta.update(
        {"driver": "GTiff", "height": out_img.shape[1], "width": out_img.shape[2], "transform": out_transform,
         "crs": from_epsg(epsg_code)})
    # pycrs.parser.from_epsg_code(epsg_code).to_proj4()
    with rasterio.open(fp, "w", **out_meta) as dest:
        dest.write(out_img)

    print("Cutting " + fp)

    return True


def dlgridfile(fn='epk10T_SHP.zip', gridurl='https://geoportaal.maaamet.ee/docs/pohikaart/epk10T_SHP.zip'):
    dlfile = requests.get(gridurl)
    if dlfile.status_code == 200:
        open(fn, 'wb').write(dlfile.content)
        return True
    return False


def dlelbtile(nr: str, size: int, tiledir: str = 'tmp') -> bool:
    # Download a tile from the Estonian Land Board server

    if not os.path.exists(tiledir):
        os.makedirs(tiledir)

    fname = tiledir + "/" + str(nr) + "_dem_" + str(size) + "m.tif"
    if not os.path.isfile(fname):
        url = "https://geoportaal.maaamet.ee/index.php?lang_id=1&plugin_act=otsing&kaardiruut=" + \
              str(nr) + "&andmetyyp=dem_" + str(size) + "m_geotiff&dl=1&f=" + str(nr) + "_dem_" + str(size) + \
              "m.tif&no_cache=5c89079809725&page_id=614"

        print("Downloading :" + url)
        dlfile = requests.get(url)
        if dlfile.status_code == 200:
            open(fname, 'wb').write(dlfile.content)
    return True


# Merge tiles together
def mergetiles(nm, flst, size, tiledir = 'tmp'):
    """

    :param nm: Filename of resulting file
    :param flst: List of tile numbers
    :param size: Resolution of the tile
    :param tiledir: Directory where the tiles go
    """
    # Check all files, download if required and merge

    if not os.path.exists(tiledir):
        os.makedirs(tiledir)

    flstr = ''
    for nr in flst:
        if len(str(nr)) > 1:
            dlelbtile(nr, size, tiledir)
            fname = tiledir + "/" + str(nr) + "_dem_" + str(size) + "m.tif"
            flstr = flstr + ' "' + fname + '"'

    cmd = 'gdal_merge.py -co COMPRESS=DEFLATE -co ZLEVEL=9 -co BIGTIFF=YES -o "' + str(nm) + '" ' + flstr
    if not os.path.isfile(str(nm)):
        subprocess.call(cmd, shell=True)
    print("Merging cluster:" + str(nm))


def dlarea(x1: float, y1: float, x2: float, y2: float, fn: str, gridshp: str = "epk10T_SHP.zip", resolution: int = 5, tiledir: str = "tmp") -> bool:
    """
    Download an area specified by coordinates

    :param x1: Min X coordinate of the area
    :param y1: Min Y coordinate of the area
    :param x2: Max X coordinate of the area
    :param y2: Max Y coordinate of the area
    :param fn: resulting file name
    :param gridshp: shapefile name of the grid
    :param resolution: resolution of the tile 1, 5 or 10
    :return: True
    """
    if resolution not in (1, 5, 10):
        print("Wrong resolution, must be etiher 1, 5, 10")
        return False

    # Open file of Estonian grid
    if not os.path.isfile(gridshp):
        print("Downloading grid file: " + gridshp)
        dlgridfile(fn = gridshp)
        return False

    if not os.path.exists(tiledir):
        os.makedirs(tiledir)

    if gridshp.lower().endswith('.zip'):
        input = fiona.open('zip://' + gridshp, layer=0)
    else:
        input = fiona.open(gridshp, 'r')

    # Find all numbers of grid cells that intersect with given box
    cells = []
    clipped = input.filter(bbox=(x1, y1, x2, y2))

    # Download all tiles
    for elem in clipped:
        cells.append(elem['properties']['NR'])

    # DL and merge all tiles
    mergetiles(fn, cells, resolution, tiledir = tiledir)

    return True


def getclosepts(pt, pts, maxrange=2000):
    """
    Get closest points to point pt from point array pts
    :param pt: Point to which the points should be closest in format {'X':X,'Y':Y}
    :param pts: List of point dicts in format {'X':X,'Y':Y}
    :param maxrange: Range of consideration
    :return: List of point dicts in format {'X':X,'Y':Y}
    """
    cpts = []

    for p in pts:
        if distance.euclidean((pt['X'], pt['Y']), (p['X'], p['Y'])) <= maxrange:
            cpts.append(p)

    for p2 in cpts:
        # Remove those points from main array
        if p2 in pts:
            pts.remove(p2)
            print("remove")
        # Get all close points to all new points
        c = getclosepts(p2, pts, maxrange)
        if type(c) is list:
            cpts = cpts + c
    return cpts


def clusterpts(pts, maxrange=2000):
    """

    :param pts: List of point dicts in format {'X':X,'Y':Y}
    :param maxrange: distance between objects to be considered one cluster
    :return:
    """
    C = []
    pts2 = copy.copy(pts)
    for pt in pts2:
        if pt in pts:
            pts.remove(pt)
            grp = getclosepts(pt, pts, maxrange)
            if type(grp) is list:
                grp.append(pt)
                C.append(grp)
                for p in grp:
                    if p in pts:
                        pts.remove(p)
            else:
                C.append([pt])
            print(len(pts))
    return C


def clt2coord(clsts, dist=2000):
    """
    Create a list of coordinates from list of point clusters

    :param clsts: List of point clusters (lists), every point is in format {'X':X,'Y':Y}
    :param dist: Radius around every point that needs to be included in the raster image
    :return: dict of coordinates defining the box
    """
    coords = []
    for c in clsts:
        xs = [x['X'] for x in c]
        ys = [x['Y'] for x in c]
        coords.append({'X1': min(xs) - dist, 'Y1': min(ys) - dist, 'X2': max(xs) + dist, 'Y2': max(ys) + dist})
    return (coords)


# Data a pandas dataframe containing X, Y fields
# Maxrange - distance between objects to be considered one cluster

def getptsmaps(data, conf):
    """

    :param data: pandas dataframe of points including coordinates with column names X, Y
    :param conf: dict of configuration variables including maxrange, resolution, gridshp
    :return: bool
    """
    if "resolution" not in conf:
        conf["resolution"] = 5

    if "gridshp" not in conf:
        conf["gridshp"] = ''

    if "tiledir" not in conf:
        conf["tiledir"] = 'tmp'

    if conf["resolution"] not in (1, 5, 10):
        print("Wrong resolution, must be either 1, 5 or 10")
        return False

    ptarr = data[["X", "Y"]].to_dict(orient="record")

    # Cluster points by distance getting n close pointsets
    if conf['maxrange']:
        areas = clusterpts(ptarr, conf['maxrange'])
        print(len([item for sublist in areas for item in sublist]))
    else:
        areas = clusterpts(ptarr)

    crd = clt2coord(areas, dist=2000)

    # Download
    n = 0
    for a in crd:
        n = n + 1
        nx1 = str(int(math.floor(a['X1'] / 1000)))
        ny1 = str(int(math.floor(a['Y1'] / 1000)))
        nx2 = str(int(math.floor(a['X2'] / 1000)))
        ny2 = str(int(math.floor(a['Y2'] / 1000)))
        mfn = conf.rasterdir + "DEM_" + nx1 + ny1 + nx2 + ny2 + '_r' + str(conf["resolution"]) + '.tif'
        if not os.path.isfile(mfn):
            dlarea(a['X1'], a['Y1'], a['X2'], a['Y2'], mfn, conf["gridshp"], conf["resolution"], conf["tiledir"])
            cropraster(mfn, a['X1'], a['Y1'], a['X2'], a['Y2'])

    return True
