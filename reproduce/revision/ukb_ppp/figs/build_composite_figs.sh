#!/bin/bash
# Rebuild the composite figures of the revision from their panels and copy them
# to paper/figs/ (local copy only; nothing is pushed to Overleaf).
#   bash code/sh_run/build_composite_figs.sh [fig2 fig3 fig4 fig5]
# Each composite is a LaTeX file next to the scripts that make its panels:
#   fig2, fig3: simulation_code/rev1/figs_compose/fig{2,3}_rev1.tex
#   fig4, fig5: code/post_cf3/figs/fig{4,5}_rev1.tex
# A composite whose panels are missing is skipped with a message.
BASE=<PROJECT_DIR>
FIGS=${@:-fig2 fig3 fig4 fig5}
for f in $FIGS; do
  case $f in fig2|fig3) d=$BASE/simulation_code/rev1/figs_compose ;; *) d=$BASE/code/post_cf3/figs ;; esac
  [ -f $d/${f}_rev1.tex ] || { echo "$f: no template"; continue; }
  ( cd $d && timeout 600 pdflatex -interaction=nonstopmode -halt-on-error ${f}_rev1.tex > ${f}_rev1.build.log 2>&1 )
  if [ $? -eq 0 ] && [ "$(pdfinfo $d/${f}_rev1.pdf 2>/dev/null | awk '/^Pages/{print $2}')" = "1" ]; then
    cp $d/${f}_rev1.pdf $BASE/paper/figs/${f}_rev1.pdf
    echo "$f: built ($(ls -la $BASE/paper/figs/${f}_rev1.pdf | awk "{printf \"%.1f MB\", \$5/1e6}")) -> paper/figs/${f}_rev1.pdf"
    rm -f $d/${f}_rev1.aux $d/${f}_rev1.log $d/${f}_rev1.build.log
  else
    echo "$f: FAILED or not one page, see $d/${f}_rev1.build.log"; grep -m3 "^!" $d/${f}_rev1.build.log
  fi
done
