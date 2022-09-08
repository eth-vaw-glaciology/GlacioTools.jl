# This tests the Antarctica download and readers.  Note that this downloads several GBs of data!
# Thus it is not run automatically.
using GlacioTools

datadir="data/antarctica"
GlacioTools.fetch_antarctica()

bm = GlacioTools.read_bedmachine(datadir, 10)[1]
bmap = GlacioTools.read_bedmap2(datadir)
basalmelt = GlacioTools.read_basal_melt_amery(dims(bmap), datadir)
lbrq = GlacioTools.read_lebrocq_flux(datadir)
ws_lakes = GlacioTools.read_wrightsiegert_lakes(datadir)
measures_gl = GlacioTools.read_gl_measures(datadir)
hogg_lakes = GlacioTools.read_hogg_lakes(datadir)
dow = GlacioTools.read_glads_dow(datadir)
mal = GlacioTools.read_malczyk_lakes(datadir)

# ## Can be verified with plotting
# using Plots; pyplot()
# plot(bm)
# plot(bmap)
# plot(basalmelt)
# plot(lbrq)
# scatter(ws_lakes[1,:], ws_lakes[2,:])
# plot(measures[1,:], measures[2,:])
# plot(hogg_lakes[:a])
# plot(dow[1])
# plot(mal[1])
