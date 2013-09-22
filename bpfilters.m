# 1:1 impedance, 50 ohm, shunt-first, 9th order, butterworth bandpass filters
# VA3DGN

function x = denorm_bp_shunt_l (l, fo, fc, rl)
  x = (rl * fc * 1000000) / (2 * pi * ((fo * 1000000) ^ 2) * l);
endfunction

function x = denorm_bp_shunt_c (c, fo, fc, rl)
  x = c / (2 * pi * fc * 1000000 * rl);
endfunction

function x = denorm_bp_series_l (l, fo, fc, rl)
  x = (rl * l) / (2 * pi * fc * 1000000);
endfunction

function x = denorm_bp_series_c (c, fo, fc, rl)
  x = (fc * 1000000) / (2 * pi * ((fo * 1000000) ^ 2) * c * rl);
endfunction

function x = butt_coeff (k, n)
  x = 2 * sin ((((2 * k) - 1) * pi) / (2 * n));
endfunction

function x = butt_table (n, t)
  if (t == 1)  # lowpass
    for k = 1:n
      x (:, k) = butt_coeff (k, n);
    endfor
  elseif (t == 2)  # highpass
    for k = 1:n
      x (:, k) = 1 / butt_coeff (k, n);
    endfor
  endif
endfunction

rs = 50;
rl = 50;
rr = rs / rl;

n = 9;
v = butt_table (n, 1);

# band edges in MHz with 3 dB point offset estimates
a = [ #low      #high   #down  #up
     1.800,     2.000,  0.05,  0.05  #160m
     3.500,     4.000,  0.05,  0.05  #80m
     5.332,     5.405,  0.05,  0.05  #60m
     7.000,     7.300,  0.05,  0.05  #40m
    10.100,    10.150,  0.05,  0.05  #30m
    14.000,    14.350,  0.05,  0.05  #20m
    18.068,    18.168,  0.05,  0.05  #17m
    21.000,    21.450,  0.05,  0.05  #15m
    24.890,    24.990,  0.05,  0.05  #12m
    28.000,    29.690,  0.05,  0.05  #10m
    50.000,    54.000,  0.05,  0.05  #6m
   144.000,   147.990,  0.05,  0.05  #2m
   220.000,   225.000,  0.05,  0.05  #125cm
   430.025,   450.000,  0.05,  0.05  #70cm
   902.000,   928.000,  0.05,  0.05  #33cm
  1240.000,  1300.000,  0.05,  0.05  #23cm
  2300.000,  2450.000,  0.05,  0.05  #13cm
];

for i = 1:size (a, 1)
  b (i, 1) = a (i, 1);                    # lower band edge
  b (i, 2) = a (i, 2);                    # upper band edge
  b (i, 3) = a (i, 1) - a (i, 3);         # lower 3 dB corner
  b (i, 4) = a (i, 2) + a (i, 4);         # upper 3 dB corner
  b (i, 5) = sqrt (b (i, 4) * b (i, 3));  # geometric centre
  b (i, 6) = b (i, 4) - b (i, 3);         # cutoff frequency / bandwidth
endfor

for i = 1:size (b, 1)
  l (i, 1) = denorm_bp_shunt_l (v (:, 1), b (i, 5), b (i, 6), rl);
  c (i, 1) = denorm_bp_shunt_c (v (:, 1), b (i, 5), b (i, 6), rl);

  l (i, 2) = denorm_bp_series_l (v (:, 2), b (i, 5), b (i, 6), rl);
  c (i, 2) = denorm_bp_series_c (v (:, 2), b (i, 5), b (i, 6), rl);

  l (i, 3) = denorm_bp_shunt_l (v (:, 3), b (i, 5), b (i, 6), rl);
  c (i, 3) = denorm_bp_shunt_c (v (:, 3), b (i, 5), b (i, 6), rl);

  l (i, 4) = denorm_bp_series_l (v (:, 4), b (i, 5), b (i, 6), rl);
  c (i, 4) = denorm_bp_series_c (v (:, 4), b (i, 5), b (i, 6), rl);

  l (i, 5) = denorm_bp_shunt_l (v (:, 5), b (i, 5), b (i, 6), rl);
  c (i, 5) = denorm_bp_shunt_c (v (:, 5), b (i, 5), b (i, 6), rl);

  l (i, 6) = denorm_bp_series_l (v (:, 6), b (i, 5), b (i, 6), rl);
  c (i, 6) = denorm_bp_series_c (v (:, 6), b (i, 5), b (i, 6), rl);

  l (i, 7) = denorm_bp_shunt_l (v (:, 7), b (i, 5), b (i, 6), rl);
  c (i, 7) = denorm_bp_shunt_c (v (:, 7), b (i, 5), b (i, 6), rl);

  l (i, 8) = denorm_bp_series_l (v (:, 8), b (i, 5), b (i, 6), rl);
  c (i, 8) = denorm_bp_series_c (v (:, 8), b (i, 5), b (i, 6), rl);

  l (i, 9) = denorm_bp_shunt_l (v (:, 9), b (i, 5), b (i, 6), rl);
  c (i, 9) = denorm_bp_shunt_c (v (:, 9), b (i, 5), b (i, 6), rl);
endfor

b
l
c
