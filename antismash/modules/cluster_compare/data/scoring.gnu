set terminal dumb size 180, 50


max_hits = 30

set xlabel 'total hits'
set ylabel 'score'
set xrange [1:max_hits]
set yrange [0:1.2]

# bonus for contiguous genes
contig_at = 0.8
contig_value = 0.8

split = 0.75      # penalty for split sections of contiguous genes
strand = 0.9      # penalty for strand mismatch of a section

base(s) = (1/contig_value)**(1/((1-contig_at)*s))
contiguous = base(max_hits)
core(x, s) = base(s)**x
f(x, s) = core(x, s) / core(s, s)
h(x, s, segments) = x > segments ? f(x,s) * split**(segments-1) : 0
g(x, s, segments, incorrect_segments) = x > segments ? h(x, s, segments) * strand**incorrect_segments : 0


plot f(x, max_hits) title sprintf('%.2f, /%d', contiguous, max_hits) dashtype 1, \
     h(x, max_hits, 2) title sprintf('%.2f, /%d, 2 segments', contiguous, max_hits) dashtype 3, \
     g(x, max_hits, 3, 1) title sprintf('%.2f, /%d, 3 segment, 1 rev strand', contiguous, max_hits) dashtype 4, \
     h(x, max_hits, 4) title sprintf('%.2f, /%d, 4 segments', contiguous, max_hits) dashtype 5, \
     g(x, max_hits, 4, 2) title sprintf('%.2f, /%d, 4 segment, 2 rev strand', contiguous, max_hits) dashtype 6, \
     h(x, max_hits, 10) title sprintf('%.2f, /%d, 10 segment', contiguous, max_hits) dashtype 7, \
     g(x, max_hits, 10, 5) title sprintf('%.2f, /%d, 10 segment, 5 rev strand', contiguous, max_hits) dashtype 8 \

#pause -1 "any key to exit"
