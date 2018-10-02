package require pmepot
package require pbctools

pbc wrap -all -center origin


volmap density [atomselect 0 "segname PROA and resid 804"] -res 0.5 -weight mass -allframes -combine avg -mol top -checkpoint 0 -minmax {{-60 -60 -120} {60 60 120}} -o d804_box.dx
pmepot -sel [atomselect 0 "not segname PROA and resid 804"] -ewaldfactor 5.0 -grid 0.5  -frames all -updatesel yes -dxfile pme_d804_box.dx -cell  {{0 0 0} { 120 0 0} { 0 120 0} { 0 0 240 }}

volmap density [atomselect 0 "segname PROA and resid 808"] -res 0.5 -weight mass -allframes -combine avg -mol top -checkpoint 0 -minmax {{-60 -60 -120} {60 60 120}} -o d808_box.dx
pmepot -sel [atomselect 0 "not segname PROA and resid 808"] -ewaldfactor 5.0 -grid 0.5  -frames all -updatesel yes -dxfile pme_d808_box.dx -cell  {{0 0 0} { 120 0 0} { 0 120 0} { 0 0 240 }}

volmap density [atomselect 0 "segname PROA and resid 926"] -res 0.5 -weight mass -allframes -combine avg -mol top -checkpoint 0 -minmax {{-60 -60 -120} {60 60 120}} -o d926_box.dx
pmepot -sel [atomselect 0 "not segname PROA and resid 926"] -ewaldfactor 5.0 -grid 0.5  -frames all -updatesel yes -dxfile pme_d926_box.dx -cell  {{0 0 0} { 120 0 0} { 0 120 0} { 0 0 240 }}

volmap density [atomselect 0 "segname PROA and resid 327"] -res 0.5 -weight mass -allframes -combine avg -mol top -checkpoint 0 -minmax {{-60 -60 -120} {60 60 120}} -o e327_box.dx
pmepot -sel [atomselect 0 "not segname PROA and resid 327"] -ewaldfactor 5.0 -grid 0.5  -frames all -updatesel yes -dxfile pme_e327_box.dx -cell  {{0 0 0} { 120 0 0} { 0 120 0} { 0 0 240 }}

volmap density [atomselect 0 "segname PROA and resid 779"] -res 0.5 -weight mass -allframes -combine avg -mol top -checkpoint 0 -minmax {{-60 -60 -120} {60 60 120}} -o e779_box.dx
pmepot -sel [atomselect 0 "not segname PROA and resid 779"] -ewaldfactor 5.0 -grid 0.5  -frames all -updatesel yes -dxfile pme_e779_box.dx -cell  {{0 0 0} { 120 0 0} { 0 120 0} { 0 0 240 }}

