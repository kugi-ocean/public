'''水温鉛直断面分布を描く'''

import numpy as np
import matplotlib.pyplot as plt  #- 描画ライブラを短縮名 plt として読み込む

### 入力
temp=np.loadtxt('DBtemp.dat')


### 格子の設定
xe = np.arange(11.0)
#- numpy arange - [0,1,2,3,...] というように1ずつ増加する値を持つ配列を作る
#-   引数を超えない数まで作成するので、ここでは0 から 10 まで要素数11の配列
ye = np.arange(11.0)*10.0

X1, Y1 = np.meshgrid(xe,ye)
#- meshgrid - 2つの1次元配列から2次元配列を作る

### 描画
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#- 図の準備

#### カラーシェード
cf1 = ax1.contourf(X1,Y1,temp,cmap='jet')
#- contourf - x,y座標を指定して temp のコンター図を描く
#-   cmap     - colormap (パレット)
plt.colorbar(cf1)
#- カラーバーを描く (カラーバー情報をcf1に付加する？)

#### 等値線図
cs1 = ax1.contour(X1,Y1,temp,20, colors='w', linewidths = 1.0)
#- 20 - 20本の等値線。
#- colors='w' - 白色の線
#- linewidths - 線の太さ
cs1.clabel(cs1.levels[::4])
#- 等値線上に数字を示す
#- [::4] - [0:終わり:4]を省略した記法。4つ飛ばしの配列となる

#### グラフ軸、ラベル (設定しなければデフォルト値が使われる)
ax1.set_xlim(0,10)  #- x軸範囲
ax1.set_ylim(100,0)  #- y軸範囲
ax1.set_title('temperature cross section (C)')
ax1.set_xlabel('distance (km)')  #- x軸ラベル
ax1.set_ylabel('depth (m)')  #- y軸ラベル

#### ディスプレイ表示とファイル保存
fig.tight_layout()  #- fig で定義されている図の余白を調整して、ラベルが重ならないように設定
plt.show()
#fig.show()  ipython で実行する場合はこちら
fig.savefig('temperature.pdf')
#plt.close()

