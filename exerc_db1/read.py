'''write.f90で作成したデータファイルを読み込む'''
import numpy as np # numpyの読み込み

nx=2  # 配列サイズ(2x3x4)を設定する
ny=3
nz=4

f=open('data.dat','rb')             # ファイルをバイナリ読み込みでオープン
a = np.fromfile(f,'d',nx*ny*nz)     # 倍精度実数を nx*ny*nz 個 a に読み込み
b = np.fromfile(f,'i8',nx*ny*nz)    # 8バイト整数を nx*ny*nz 個 b に読み込み
f.close()                           # ファイルを閉じる

a = a.reshape((nx,ny,nz),order='F') # a を Fortranの順の3次元配列に変形
b = b.reshape((nx,ny,nz),order='F') # b を Fortranの順の3次元配列に変形

print(a) # aの中身の表示 
print(b) # bの中身の表示
