function [region] = findRegion(rec, mask)

dphiRec1 = conv2(rec, [1 -1]);      dphiRec1 = dphiRec1(:, 1:end-1) ;
dphiRec2 = conv2(rec, [1 -1]);      dphiRec2 = dphiRec2(:, 1:end-1) ;
dphiRec = (((dphiRec1 + dphiRec2) < -2.0) .* 1 + mask) .* mask;
region = scanline(dphiRec);
