function X = get_x_pbc(data)
  % trt + age + sex + edema + bili + albumin + platelet
  trt = data.trt;
  age = data.age;
  edema5 = data.edema == 0.5;
  edema1 = data.edema == 1;
  bili = data.bili;
  albumin = data.albumin;
  platelet = data.platelet;
  X = [trt, age, edema5, edema1, bili, albumin, platelet];
end
