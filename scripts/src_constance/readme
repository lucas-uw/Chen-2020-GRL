## Description of scripts


step1.identify_AR.argv.py and step1.identify_AR.argv.ext.py (for extended period) :
identify the AR patch (abs and p85). This requires a combined input from NARR and WRF (regridded to NARR grids), so different from step3


step2.adjust_NARR_ARtag_using_WRF.bp.py and step2.adjust_NARR_ARtag_using_WRF.py (for extended period) :
modify the AR pixels using the WRF IWV and uvIVT

step3.identify_NARR_AR.argv.py and step3.identify_NARR_AR.py (for extended period) :
identify the AR patch in the NARR data


step4.extract_ARstats.py
compute the intermediate data

step5.moisture_decomposition.py
compute the moisture component within the AR area (so IVT is only the west boundary and AR region)

step6.moisture_decomposition.no_AR_boundary.py
compute the moisture component over the ocean (so IVT is the whole western boundary, plus north and south boundary of the ocean)

step7.extract_IVT_SST.full_ocean.py
compute the IVT and IWV over the ocean, for analyzing the sensitivity of IWV/IVT to SST overa ocean (rather than in step4 where the landfalling part of IWV/IVT is analyzed)
