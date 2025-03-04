## Zhilong & CMikolaitis @ USA

from pylib import *
import pandas as pd

## Config----------------------------------------------------------------------
source = 'mesh82_utm18n'    #Mesh name not extension
shapefile = False           #Output a shapefile?

## Core functions--------------------------------------------------------------
# Convert 2d mesh to gr3 then open
inname  = source + '.2dm'
outname = source + '.gr3'
sms2grd(inname,outname)

gd_org = read_schism_hgrid(outname)
gd_org.dp = -gd_org.dp
gd_org.write_hgrid(outname)

gd_org.compute_all()

# Calculate CFLs & put into df
t = 100  # seconds
i = 0
CFLs = []; hs = []; xs = []
for i in range(gd_org.ne):
    A = gd_org.area[i]      # area of element
    x = sqrt(A*4/pi)        # meters
    ni = gd_org.elnode[i]   # nodes at element i
    if gd_org.i34[i]==3:
        h = mean(gd_org.dp[ni[0:3]])
    else:
        h = mean(gd_org.dp[ni])
    if h >= 0.1:
       CFL = (sqrt(9.81*h)+1)*t/x   # calculate CFL number h >= 0.1 m
    else:
        CFL = -999                   # calculate CFL number  h < 0.1 m
    CFLs.append(CFL)
    hs.append(h)
    xs.append(x)
df = pd.DataFrame({'CFL':CFLs})
test = pd.DataFrame({'h':hs,'x':xs,'CFL':CFLs})

# Create masks
lCFL,mCFL,uCFL,sCFL = 0.5,1.0,5.0,-999
shal = df[df==sCFL].dropna().index.values.tolist()
low  = df[(df<=lCFL) & (df!=sCFL)].dropna().index.values.tolist()
mid  = df[(df<=mCFL) & (df>lCFL)].dropna().index.values.tolist()
hi   = df[df>uCFL].dropna().index.values.tolist()

# Mask
shal  = pd.DataFrame({'x':gd_org.xctr[shal],'y':gd_org.yctr[shal],'CFL':test['CFL'].iloc[shal]})
lowDF = pd.DataFrame({'x':gd_org.xctr[low],'y':gd_org.yctr[low],'CFL':test['CFL'].iloc[low]})
midDF = pd.DataFrame({'x':gd_org.xctr[mid],'y':gd_org.yctr[mid],'CFL':test['CFL'].iloc[mid]})
hiDF  = pd.DataFrame({'x':gd_org.xctr[hi],'y':gd_org.yctr[hi],'CFL':test['CFL'].iloc[hi]})

# Get stats
print(test.describe())

# Get plots
msize = 4
figure(figsize=[80,40])
gd_org.plot(fmt=0,ec='k',lw=0.05,clim=[0,100],ticks=11)
plot(shal.x,shal.y,'g.',ms=msize)
plot(lowDF.x,lowDF.y,'b.',ms=msize)
plot(midDF.x,midDF.y,'y.',ms=msize)
plot(hiDF.x,hiDF.y,'r.',ms=msize)
title("CFLs: t=%i, Blue < %0.1f, Yellow < %0.1f, Red > %0.1f, Green: Too Shallow" %(t,lCFL,mCFL,uCFL))
savefig('meshtest.png'); close()

# Make shapefile
def shpmkr(source,shpname):
    C           = zdata()
    C.type      = 'Point'
    C.prj       = 'epsg:26918'
    C.xy        = c_[source.x,source.y]
    C.attname   = ['CFL']
    C.attvalue  = source.CFL
    write_shp(shpname,C)
    
# Get shapefile
if shapefile:
    shpmkr(lowDF,low_CFLs)