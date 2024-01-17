1. Check the packages_list.md to see if all the necessary packages are installed on your system.
2. The steps of sh file excution: QC > BWA > markduplicate > bqsr > summary
3. For sample only runs on a single lane, mapping and markduplicate would be a bit different from general example files. Please use: bwa_mkdp_single_lane_example.sh
4. For sample runs on multiple lanes, raw reads from each lane have to be conducted separately from QC to bwa. The files would be merged at markduplicate.
