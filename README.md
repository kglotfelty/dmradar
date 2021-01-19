# `dmradar`

This is the polar equivalent of the CIAO 
[`dmnautilus`](https://cxc.cfa.harvard.edu/ciao/ahelp/dmnautilus.html) 
tool.

It works by subdividing the image into 4 pie-shaped wedges and checking
if the SNR in each wedge meets the threshold criteria.  If so, then each wedge
is divided in 4 (radii by 2, angle divided by 2) and the process is repeated.

By default, the wedges are circular: `pie` shapes.  Users can also
use elliptical wedges, `epanda` , or box wedges, `bpanda`.  There is also
a `box` mode which emulates the behavior of the original `dmnautilus` tool.


> **Where does `*panda` come from?**
> 
> `panda` is from `SAOImage ds9`.  It means: _pie and annulus_.

