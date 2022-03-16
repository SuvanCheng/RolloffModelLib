# RolloffModelLib



# RolloffModel

## 运行环境
> win11操作系统下，vs2022、nodejs、IDEA

## 前言
本系统实现通过碾压车碾压数据拟合大坝建设动态模型，总体实现三个功能：建模、转模型格式、写czml文件。

## 流程图

①	总流程图

![image](https://user-images.githubusercontent.com/24319035/158499989-7b803996-0d14-4181-a47e-d36ca7ab7392.png)

② 建模流程图

![image](https://user-images.githubusercontent.com/24319035/158500002-327731c1-59b5-4932-bdbc-36d61e76c7e6.png)

## 核心原理：凸包(Convex Hull)构造算法Graham扫描法
① 找到所有点P0,1,...,N−1P0,1,...,N−1中最下方的点，记为PLPL；
② 计算所有其他的点Pi(i≠L)Pi(i≠L) 与 PLPL构成的向量PLPi−→−−PLPi→相对于水平轴的夹角。因为所有的点都在该PLPL上方，因此向量的取值范围为(0,180)(0,180) ，所以可以用余切值代替角度值；
③ 对所有其他的点按照第2步算出的角度进行排序，且PLPL为排序后的数组的第0位；
④ 从点PLPL开始，依此连接每一个点（已经排序过），每连接一个点检测连线的走向是否是逆时针的，如果是则留下该点的前一个点，反之去除前一个点，使之与前面第二个点直接连接，继续这一检测，直到是逆时针或者所有点都被检测过为止。


## 使用说明
本程序提供cpp源码和Java程序调用接口，总共两个程序使用方式大体一样
```java
//SingleModel
public native void ModelC(
    double[] x,double[] y,double [] z, //表示数据x轴坐标y轴坐标z轴坐标
    //double[] time, //时间戳，在程序中与xyz轴坐标用数组存储
    String path, //程序输出的文件夹，必须存在，用字符串存储
    double dz, //z轴精度，程序按照z轴精度依次次探查
    int lim, //输入一单位z轴值的平面上允许的最少个数，例如10，如果这一单位z轴值的平面上有的点个数小于输入值，则删除这一平面上的点
    int add_float, //数值是判断相连几个面组成的层内是否点云个数过少
    int min_point_num, //输入一单位z轴值的平面上允许的最少个数，这个值用于判断相连几个面组成的层内是否点云个数过少
    double max_face, //一个值用来判断这一层的点云是否距离过大
    double min_face, //一个值用来判断这一层重建的是否为道路
    //String longitude,String latitude,String heightDeo, //经度、纬度和高程，用来写czml文件
    //int timegap, //设定用几天的时间内的数据一起建模
    double mov_x,double mov_y //两个值用来将模型移动到坐标中心点。如果输入一个为0，则默认都使用模型中值进行平移
    );

//TimeModel
public native void ModelC(
    double[] x,double[] y,double [] z,double[] time, //表示数据x轴坐标y轴坐标z轴坐标和时间戳，在程序中用数组存储。
    String path, //程序输出的文件夹，必须存在，用字符串存储
    double dz, //z轴精度，程序按照z轴精度依次次探查
    int lim, //输入一单位z轴值的平面上允许的最少个数，例如10，如果这一单位z轴值的平面上有的点个数小于输入值，则删除这一平面上的点
    int add_float, //数值是判断相连几个面组成的层内是否点云个数过少
    int min_point_num, //输入一单位z轴值的平面上允许的最少个数，这个值用于判断相连几个面组成的层内是否点云个数过少
    double max_face, //一个值用来判断这一层的点云是否距离过大
    double min_face, //一个值用来判断这一层重建的是否为道路
    String longitude,String latitude,String heightDeo, //经度、纬度和高程，用来写czml文件
    int timegap, //设定用几天的时间内的数据一起建模
    double mov_x,double mov_y //两个值用来将模型移动到坐标中心点。如果输入一个为0，则默认都使用模型中值进行平移
    );
```
通过在java程序中调用该函数，即可运行


