from astroquery.vsa import Vsa
from astropy.coordinates import Angle, SkyCoord
import astropy.units as u
import pandas as pd
from io import StringIO
import re
import requests

#%%
Vsa.ROW_LIMIT = -1


columns = [
    'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
    'pSaturated', 'yPetroMag', 'yPetroMagErr', 'yPsfMag', 'yPsfMagErr', 'ySerMag2D', 'ySerMag2DErr', 'yAperMag3',
    'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yAperMag6', 'yAperMag6Err', 'yAperMagNoAperCorr3',
    'yAperMagNoAperCorr4', 'yAperMagNoAperCorr6', 'yEll', 'yPA', 'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
    'jPetroMag', 'jPetroMagErr', 'jPsfMag', 'jPsfMagErr', 'jSerMag2D', 'jSerMag2DErr', 'jAperMag3',
    'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jAperMag6', 'jAperMag6Err', 'jAperMagNoAperCorr3',
    'jAperMagNoAperCorr4', 'jAperMagNoAperCorr6', 'jEll', 'jPA', 'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
    'hPetroMag', 'hPetroMagErr', 'hPsfMag', 'hPsfMagErr', 'hSerMag2D', 'hSerMag2DErr', 'hAperMag3',
    'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hAperMag6', 'hAperMag6Err', 'hAperMagNoAperCorr3',
    'hAperMagNoAperCorr4', 'hAperMagNoAperCorr6', 'hEll', 'hPA', 'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',
    'ksPetroMag', 'ksPetroMagErr', 'ksPsfMag', 'ksPsfMagErr', 'ksSerMag2D', 'ksSerMag2DErr', 'ksAperMag3',
    'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksAperMag6', 'ksAperMag6Err', 'ksAperMagNoAperCorr3',
    'ksAperMagNoAperCorr4', 'ksAperMagNoAperCorr6', 'ksEll', 'ksPA', 'ksErrBits', 'ksClass', 'ksClassStat', 'ksppErrBits'
]


coord = SkyCoord(ra=0*u.deg, dec=-90*u.deg)

table2 = Vsa.query_region(
    coord,
    radius=2*u.deg,
    programme_id='VHS',
    database='VHSDR6',
    attributes=columns,
    constraints="priOrSec=0"
)
#%%
print(len(table2))

#%%
columns = [
    'sourceID', 'ra', 'dec', 'l', 'b', 'priOrSec', 'mergedClassStat', 'mergedClass', 'pStar', 'pGalaxy', 'pNoise',
    'pSaturated', 'yPetroMag', 'yPetroMagErr', 'yPsfMag', 'yPsfMagErr', 'ySerMag2D', 'ySerMag2DErr', 'yAperMag3',
    'yAperMag3Err', 'yAperMag4', 'yAperMag4Err', 'yAperMag6', 'yAperMag6Err', 'yAperMagNoAperCorr3',
    'yAperMagNoAperCorr4', 'yAperMagNoAperCorr6', 'yEll', 'yPA', 'yErrBits', 'yClass', 'yClassStat', 'yppErrBits',
    'jPetroMag', 'jPetroMagErr', 'jPsfMag', 'jPsfMagErr', 'jSerMag2D', 'jSerMag2DErr', 'jAperMag3',
    'jAperMag3Err', 'jAperMag4', 'jAperMag4Err', 'jAperMag6', 'jAperMag6Err', 'jAperMagNoAperCorr3',
    'jAperMagNoAperCorr4', 'jAperMagNoAperCorr6', 'jEll', 'jPA', 'jErrBits', 'jClass', 'jClassStat', 'jppErrBits',
    'hPetroMag', 'hPetroMagErr', 'hPsfMag', 'hPsfMagErr', 'hSerMag2D', 'hSerMag2DErr', 'hAperMag3',
    'hAperMag3Err', 'hAperMag4', 'hAperMag4Err', 'hAperMag6', 'hAperMag6Err', 'hAperMagNoAperCorr3',
    'hAperMagNoAperCorr4', 'hAperMagNoAperCorr6', 'hEll', 'hPA', 'hErrBits', 'hClass', 'hClassStat', 'hppErrBits',
    'ksPetroMag', 'ksPetroMagErr', 'ksPsfMag', 'ksPsfMagErr', 'ksSerMag2D', 'ksSerMag2DErr', 'ksAperMag3',
    'ksAperMag3Err', 'ksAperMag4', 'ksAperMag4Err', 'ksAperMag6', 'ksAperMag6Err', 'ksAperMagNoAperCorr3',
    'ksAperMagNoAperCorr4', 'ksAperMagNoAperCorr6', 'ksEll', 'ksPA', 'ksErrBits', 'ksClass', 'ksClassStat', 'ksppErrBits'
]

joined_cols = ','.join(columns)

BASE_URL = "http://vsa.roe.ac.uk:8080/vdfs/WSASQL"

params = {
    "programmeID": 110,
    "database": "VHSDR7",
    "archive": "VSA",
    "formaction": "region",
    "from": "source",
    "ra": 0,
    "dec": -90,
    "sys": "J",
    "name": "",
    "radius": 90,
    "xSize": "",
    "ySize": "",
    "boxAlignment": "RADec",
    "emailAddress": "",
    "format": "CSV",
    "compress": "",
    "rows": 30,
    "select": joined_cols,
    "where": "priOrSec=0",
}

response = requests.get(BASE_URL, params=params)
response.raise_for_status()

match = re.search(r'href="(http://vsa\.roe\.ac\.uk/tmp/tmp_sql/results[^"]+\.csv)"', response.text)
if not match:
    raise ValueError("Could not find CSV download link in response")

csv_url = match.group(1)
print(f"Found CSV URL: {csv_url}")

csv_response = requests.get(csv_url)
csv_response.raise_for_status()

df = pd.read_csv(StringIO(csv_response.text), comment="#")
print(df.shape)

#%%
print(df[:10])