#include <iostream>

#include "SingleModel.h"
#include <iostream>
#include <stdio.h>

#ifdef WIN32
#include <io.h>
#else
#if (defined __ARM__) || (defined __APPLE__)
#include <sys/uio.h>
#else
#include <sys/io.h>
#endif
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstring>
#include <string>
#include <chrono>
#include <cstdio>
#include <math.h>
#pragma warning(disable:4996)
using namespace std;//double比较精度问题
using std::stringstream;
#define EPS 1e-6

//点结构
class Point
{
public:
    double x, y, z, time;
    int num;
    Point(double x = 0, double y = 0, double z = 0, double time = 0, int num = 0) :x(x), y(y), z(z), time(time), num(num) {}
    Point operator + (Point a)
    {
        return Point(a.x + x, a.y + y);
    }
    Point operator - (Point a)
    {
        return Point(x - a.x, y - a.y);
    }
    bool operator < (const Point& a) const
    {
        if (x == a.x)
            return y < a.y;
        return x < a.x;
    }
};

//面结构（三点构成面）
struct Tri
{
    int a;
    int b;
    int c;
};

//时间和模型绑定结构体
struct num_day
{
    int num_dd;
    double day;
};

//实层结构（记录相连实面）
struct floatz {
    double under;
    double up;
};

typedef Point Vector;

//计算凸包前置函数
double cross(Vector a, Vector b)
{
    return a.x * b.y - a.y * b.x;
}
double dot(Vector a, Vector b)
{
    return a.x * b.x + a.y * b.y;
}

//判断连线是否组成凸包
bool isclock(Point p0, Point p1, Point p2)
{
    Vector a = p1 - p0;
    Vector b = p2 - p0;
    if (cross(a, b) < -EPS) return true;
    return false;
}
double getDistance(Point a, Point b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}
typedef vector<Point> Polygo;

//计算凸包点
Polygo andrewScan(Polygo s)
{
    Polygo u, l;
    if (s.size() < 3) return s;
    sort(s.begin(), s.end());
    u.push_back(s[0]);
    u.push_back(s[1]);
    l.push_back(s[s.size() - 1]);
    l.push_back(s[s.size() - 2]);
    //printf("l[n-2]:%.2f %.2f\nl[n-1]:%.2f %.2f\n", l[l.size() - 2].x, l[l.size() - 2].y, l[l.size() - 1].x, l[l.size() - 1].y);
    for (int i = 2; i < s.size(); i++)
    {
        for (int n = u.size(); n >= 2 && isclock(u[n - (int)2], u[n - (int)1], s[i]) != true; n--)
        {
            //cout << u[n - 2].x << ' ' << u[n - 2].y << '\n' << u[n - 1].x << u[n - 1].y << endl;
            u.pop_back();
        }
        u.push_back(s[i]);
    }
    for (int i = s.size() - 3; i >= 0; i--)
    {
        //cout << i << endl;
        for (int n = l.size(); n >= 2 && isclock(l[n - (int)2], l[n - (int)1], s[i]) != true; n--)
        {
            //cout << i << endl;
            // printf("del:\nl[n-2]:%.2f %.2f\nl[n-1]:%.2f %.2f\n", l[n - 2].x, l[n - 2].y, l[n - 1].x, l[n - 1].y);

            l.pop_back();
        }

        l.push_back(s[i]);
    }
    //    for(auto &p : u) printf("%.2f %.2f\n",p.x,p.y);
    //    printf("yes\n");
    //    for(auto &p : l) printf("%.2f %.2f\n",p.x,p.y);

    for (int i = 1; i < u.size() - 1; i++) l.push_back(u[i]);
    return l;
}
//排序算法
bool  zsortup(const  Point& s1, const  Point& s2)
{
    return  s1.z < s2.z;
}
bool  zsortdown(const  Point& s1, const  Point& s2)
{
    return  s1.z > s2.z;
}
bool  xsortup(const  Point& s1, const  Point& s2)
{
    return  s1.x < s2.x;
}
bool  xsortdown(const  Point& s1, const  Point& s2)
{
    return  s1.x > s2.x;
}
bool  ysortup(const  Point& s1, const  Point& s2)
{
    return  s1.y < s2.y;
}
bool  ysortdown(const  Point& s1, const  Point& s2)
{
    return  s1.y > s2.y;
}
bool  tsortup(const  Point& s1, const  Point& s2)
{
    return  s1.time < s2.time;
}
bool  fsortup(const  floatz& s1, const  floatz& s2)
{
    return  s1.under < s2.under;
}
//清空少点面
void clearempty(Polygo& fileroot, vector<double>& kongdown, int& time_b, int& time_r, int lim, Polygo& file2, double down, double dz)
{
    Polygo filein;
    filein.assign(fileroot.begin() + time_b, fileroot.begin() + time_r);
    sort(filein.begin(), filein.end(), zsortup);
    int read_num = 0;//filein读取位置
    while (read_num < filein.size())
    {
        Polygo zmid;
        while (filein[read_num].z + EPS < down + dz)
        {
            zmid.push_back(filein[read_num]);
            read_num++;
            if (read_num >= filein.size())
                break;
        }
        if (zmid.size() < lim)
        {
            kongdown.push_back(down);
            down += dz;
            zmid.clear();
            continue;
        }
        down += dz;
        file2.insert(file2.end(), zmid.begin(), zmid.end());
        zmid.clear();
    }
    filein.clear();
}
//合并计算层
void creatfull(Polygo& file2, vector<double>& kongdown, vector<floatz>& full, double& dz)
{
    if (file2.front().z < kongdown.front() && fabs(file2.front().z + dz - kongdown.front()) >= EPS)
    {
        floatz mid;
        mid.under = file2.front().z;
        mid.up = kongdown.front() - dz;
        full.push_back(mid);
    }
    for (int i = 0; i < kongdown.size() - 1; i++)
    {
        if (fabs(kongdown[i] + dz - kongdown[i + 1]) >= EPS)
        {
            floatz mid;
            mid.under = kongdown[i] + dz;
            mid.up = kongdown[i + 1] - dz;
            full.push_back(mid);
        }
    }
    if (kongdown.back() < file2.back().z && fabs(kongdown.back() + dz - file2.back().z) >= EPS)
    {
        floatz mid;
        mid.under = kongdown.back() + dz;
        mid.up = file2.back().z;
        full.push_back(mid);
    }

    vector<int> Toosmall;
    for (int i = 0; i < full.size() - 1; i++)
    {
        if (fabs(full[i].up - full[i].under) <= dz + EPS)
        {
            if (fabs(full[i + 1].under - full[i].up) <= 2 * dz + EPS)
            {
                Toosmall.push_back(i);
            }
        }
    }
    if (!Toosmall.empty()) {
        for (int i = 0; i < Toosmall.size() - 1; i++)
        {
            int mid = i;
            while (Toosmall[mid] + 1 == Toosmall[mid + 1])
            {
                mid++;
                if (mid == Toosmall.size() - 1)
                    break;
            }
            floatz midd;
            midd.under = full[Toosmall[i]].under;
            midd.up = full[Toosmall[mid] + 1].up;
            full.push_back(midd);
            i = mid;
        }
        if (Toosmall.size() == 1)
        {
            floatz midd;
            midd.under = full[Toosmall.back()].under;
            midd.up = full[Toosmall.back() + 1].up;
            full.push_back(midd);
        }
    }
    if (!Toosmall.empty()) {
        int del_n = 0;
        for (int i = 0; i < Toosmall.size() - 1; i++)
        {
            int lll = Toosmall[i] - del_n;
            full.erase(full.begin() + lll);
            del_n++;
            if (Toosmall[i] + 1 != Toosmall[i + 1])
            {
                int kkk = Toosmall[i] + 1 - del_n;
                full.erase(full.begin() + kkk);
                del_n++;
            }
        }
        int lll = Toosmall.back() - del_n;
        full.erase(full.begin() + lll);
        del_n++;
        int kkk = Toosmall.back() + 1 - del_n;
        full.erase(full.begin() + kkk);
        del_n++;
    }
}
//求中间底面
double daymaxz(Polygo& file2, double& day_min_z, double& day_max_z, double& down)
{
    double day_mid_max_z = file2.front().z;
    for (int i = 1; i < file2.size(); i++)
    {
        if (day_min_z > file2[i].z)
        {
            day_min_z = file2[i].z;
        }
        if (day_mid_max_z < file2[i].z)
        {
            day_mid_max_z = file2[i].z;
        }
    }
    if (day_min_z > day_max_z)
    {
        down = day_max_z;
    }
    else
    {
        down = day_min_z;
    }
    return day_mid_max_z;
}
// 写obj文件
void writeOBJ(Polygo& Out, ofstream& myout, int& num, double& dz, double& down)
{
    string k = " ";
    string v = "v";
    string f = "f";
    string h = "\n";
    Out.pop_back();
    int num_vn = 1;
    Polygo Out2;
    Out2.assign(Out.begin(), Out.end());
    if (fabs(Out[0].z - down) > 10)
    {
        down = Out[0].z - dz;
    }
    for (int i = 0; i < Out.size(); i++)
    {
        myout << setiosflags(ios::fixed) << v + k << Out[i].x << k << Out[i].z << k << Out[i].y << h;
        Out[i].num = num;
        num++;
    }
    for (int i = 0; i < Out2.size(); i++)
    {
        Out2[i].z = down;
        myout << setiosflags(ios::fixed) << v + k << Out2[i].x << k << down << k << Out2[i].y << h;
        Out2[i].num = num;
        num++;
    }
    myout << setiosflags(ios::fixed) << "vn 0 1 0" << "\n";
    num_vn++;
    myout << setiosflags(ios::fixed) << "vn 0 -1 0" << "\n";
    num_vn++;
    for (int i = 1; i < Out.size() - 1; i++)
    {
        double c = (Out[i].x - Out[0].x) * (Out[i + 1].y - Out[0].y) - (Out[i + 1].x - Out[0].x) * (Out[i].y - Out[0].y);
        if (c >= 0)
        {
            myout << setiosflags(ios::fixed) << f << k << Out[0].num << "//" << num_vn - 2 << k << Out[i].num << "//" << num_vn - 2 << k << Out[i + 1].num << "//" << num_vn - 2 << h;
            myout << setiosflags(ios::fixed) << f << k << Out[0].num << "//" << num_vn - 1 << k << Out[i + 1].num << "//" << num_vn - 1 << k << Out[i].num << "//" << num_vn - 1 << h;
        }
        else
        {
            myout << setiosflags(ios::fixed) << f << k << Out[0].num << "//" << num_vn - 2 << k << Out[i + 1].num << "//" << num_vn - 2 << k << Out[i].num << "//" << num_vn - 2 << h;
            myout << setiosflags(ios::fixed) << f << k << Out[0].num << "//" << num_vn - 1 << k << Out[i].num << "//" << num_vn - 1 << k << Out[i + 1].num << "//" << num_vn - 1 << h;
        }
    }
    myout << f;
    for (int i = 1; i < Out2.size() - 1; i++)
    {
        double c = (Out2[i].x - Out2[0].x) * (Out2[i + 1].y - Out2[0].y) - (Out2[i + 1].x - Out2[0].x) * (Out2[i].y - Out2[0].y);
        if (c >= 0)
        {
            myout << setiosflags(ios::fixed) << f << k << Out2[0].num << "//" << num_vn - 2 << k << Out2[i + 1].num << "//" << num_vn - 2 << k << Out2[i].num << "//" << num_vn - 2 << h;
            myout << setiosflags(ios::fixed) << f << k << Out2[0].num << "//" << num_vn - 1 << k << Out2[i].num << "//" << num_vn - 1 << k << Out2[i + 1].num << "//" << num_vn - 1 << h;
        }
        else
        {
            myout << setiosflags(ios::fixed) << f << k << Out2[0].num << "//" << num_vn - 2 << k << Out2[i].num << "//" << num_vn - 2 << k << Out2[i + 1].num << "//" << num_vn - 2 << h;
            myout << setiosflags(ios::fixed) << f << k << Out2[0].num << "//" << num_vn - 1 << k << Out2[i + 1].num << "//" << num_vn - 1 << k << Out2[i].num << "//" << num_vn - 1 << h;
        }
    }
    for (int i = 0; i < Out.size() - 1; i++)
    {
        double a = (Out2[i + 1].y - Out2[i].y) * (Out[i + 1].z - Out2[i].z) - (Out[i + 1].y - Out2[i].y) * (Out2[i + 1].z - Out2[i].z);
        double b = (Out2[i + 1].z - Out2[i].z) * (Out[i + 1].x - Out2[i].x) - (Out[i + 1].z - Out2[i].z) * (Out2[i + 1].x - Out2[i].x);
        double c = (Out2[i + 1].x - Out2[i].x) * (Out[i + 1].y - Out2[i].y) - (Out[i + 1].x - Out2[i].x) * (Out2[i + 1].y - Out2[i].y);
        myout << setiosflags(ios::fixed) << "vn" << k << a << k << b << k << c << h;
        num_vn++;
        myout << f + k << Out2[i].num << "//" << num_vn - 1 << k << Out2[i + 1].num << "//" << num_vn - 1 << k << Out[i + 1].num << "//" << num_vn - 1 << k << Out[i].num << "//" << num_vn - 1 << h;
        myout << setiosflags(ios::fixed) << "vn" << k << -a << k << -b << k << -c << h;
        num_vn++;
        myout << f + k << Out[i].num << "//" << num_vn - 1 << k << Out[i + 1].num << "//" << num_vn - 1 << k << Out2[i + 1].num << "//" << num_vn - 1 << k << Out2[i].num << "//" << num_vn - 1 << h;
    }
    double a = (Out2.front().y - Out2.back().y) * (Out.front().z - Out2.back().z) - (Out.front().y - Out2.back().y) * (Out2.front().z - Out2.back().z);
    double b = (Out2.front().z - Out2.back().z) * (Out.front().x - Out2.back().x) - (Out.front().z - Out2.back().z) * (Out2.front().x - Out2.back().x);
    double c = (Out2.front().x - Out2.back().x) * (Out.front().y - Out2.back().y) - (Out.front().x - Out2.back().x) * (Out2.front().y - Out2.back().y);
    myout << setiosflags(ios::fixed) << "vn" << k << -a << k << -b << k << -c << h;
    num_vn++;
    myout << f + k << Out.back().num << "//" << num_vn - 1 << k << Out.front().num << "//" << num_vn - 1 << k << Out2.front().num << "//" << num_vn - 1 << k << Out2.back().num << "//" << num_vn - 1 << h << endl;
    myout << setiosflags(ios::fixed) << "vn" << k << a << k << b << k << c << h;
    num_vn++;
    myout << f + k << Out2.back().num << "//" << num_vn - 1 << k << Out2.front().num << "//" << num_vn - 1 << k << Out.front().num << "//" << num_vn - 1 << k << Out.back().num << "//" << num_vn - 1 << h << endl;
}
//转格式
void obj2gltf(string path, int num_time_name)
{
    system("npm install -g obj2gltf");
    string sys1 = "obj2gltf -i ";
    string sys2 = " -o ";
    for (int i = 1; i < num_time_name; i++)
    {
        string sys = sys1 + path + to_string(i) + ".obj" + sys2 + path + to_string(i) + ".gltf";
        system(sys.c_str());
    }
}
//时间戳转时间
std::tm* gettm(long long timestamp)
{
    auto milli = timestamp + (long long)8 * 60 * 60 * 1000; //此处转化为东八区北京时间，如果是其它时区需要按需求修改
    auto mTime = std::chrono::milliseconds(milli);
    auto tp = std::chrono::time_point<std::chrono::system_clock, std::chrono::milliseconds>(mTime);
    auto tt = std::chrono::system_clock::to_time_t(tp);
    std::tm* now = std::gmtime(&tt);
    printf("%4d年%02d月%02d日 %02d:%02d:%02d\n", now->tm_year + 1900, now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
    return now;
    // return NULL;
}
//时间规范字符串
string DayString(std::tm* day)
{
    // 准备根据格式造字符串流
    stringstream fmt;   /* 或者使用 ostringstream */
    // 造字符串流
    fmt << day->tm_year + 1900 << "-" << setw(2) << setfill('0') << day->tm_mon + 1 << "-" << setw(2) << day->tm_mday << "T" << setw(2) << day->tm_hour << ":" << setw(2) << day->tm_min << ":" << setw(2) << day->tm_sec << ".0z";
    string targetString = fmt.str();
    return targetString;
}

JNIEXPORT void JNICALL Java_SingleModel_ModelC
        (JNIEnv* env, jobject ob, jdoubleArray jx, jdoubleArray jy, jdoubleArray jz, jstring jpath, jdouble jdz, jint jlim, jint jadd_float, jint jmin_point_num, jdouble jmax_face, jdouble jmin_face,  jdouble mov_x, jdouble mov_y) {
    int num_time_name = 1;
    double day_min_z = 0;
    double day_max_z = 0;
    vector<num_day> daylist;//时间与编码
    string file_path = "/Users/suvancheng/Desktop/BIDR/JNI";
    string path = "/Users/suvancheng/Desktop/BIDR/JNI";//输出位置，输出准备
    const char* jcpath;
    jcpath = env->GetStringUTFChars(jpath, JNI_FALSE);
    path = jcpath;
    string bu = "/"; //注意：Windows写"\\"，macOS写"/"
    double dz = 0.001;//切割高程
    dz = jdz;
    int lim = 10;//高程内最少的点数
    lim = jlim;
    double geodown = (double)0;
    int add_float = 2;
    add_float = jadd_float;
    int min_point_num = 10;
    min_point_num = jmin_point_num;
    double max_face = 150;
    max_face = jmax_face;
    double min_face = 50;
    int timegap = 1;
    min_face = jmin_face;
    string longitude = "116.6445206464275";
// const char* jclongitude;
//jclongitude = env->GetStringUTFChars(jlongitude, false);
//longitude = jclongitude;
    string latitude = "33.24206160289476";
//const char* jclatitude;
//jclatitude = env->GetStringUTFChars(jlatitude, false);
//latitude = jclatitude;
    string heightDeo = "19.5";
//const char* jcheightDeo;
//jcheightDeo = env->GetStringUTFChars(jheightDeo, false);
//heightDeo = jcheightDeo;
    Polygo fileroot;//原始点云
    double min_x = 0;
    double min_y = 0;
    double min_z = 0;
    double min_t = 0;
    double max_x = 0;
    double max_y = 0;
    double max_z = 0;
    double max_t = 0;
//------------读取点云------------------
    double* x = env->GetDoubleArrayElements(jx, 0);
    double* y = env->GetDoubleArrayElements(jy, 0);
    double* z = env->GetDoubleArrayElements(jz, 0);
//double* time = env->GetDoubleArrayElements(jtime, 0);
    int num_jx = env->GetArrayLength(jx);
    int num_jy = env->GetArrayLength(jy);
    int num_jz = env->GetArrayLength(jz);
//int num_jtime = env->GetArrayLength(jtime);
    int jmin = num_jx;
    if (jmin > num_jy)
        jmin = num_jy;
    if (jmin > num_jz)
        jmin = num_jz;
//if (jmin > num_jtime)
//    jmin = num_jtime;
    for (int i = 0; i < jmin; i++)
    {
        Point mid;
        mid.x = x[i];
        mid.y = y[i];
        mid.z = z[i];
        mid.time = 0;
        mid.num = 0;
        fileroot.push_back(mid);
    }
    min_x = fileroot[0].x;
    min_y = fileroot[0].y;
    min_z = fileroot[0].z;
    min_t = fileroot[0].time;
    max_x = fileroot[0].x;
    max_y = fileroot[0].y;
    max_z = fileroot[0].z;
    max_t = fileroot[0].time;
    for (int i = 1; i < fileroot.size(); i++)//找最小值
    {
        if (min_x > fileroot[i].x)
            min_x = fileroot[i].x;
        if (max_x < fileroot[i].x)
            max_x = fileroot[i].x;
        if (min_y > fileroot[i].y)
            min_y = fileroot[i].y;
        if (max_y < fileroot[i].y)
            max_y = fileroot[i].y;
        if (min_z > fileroot[i].z)
            min_z = fileroot[i].z;
        if (max_z < fileroot[i].z)
            max_z = fileroot[i].z;
    }
    double mid_x = (float)(max_x + min_x) / 2.0;
    double mid_y = (float)(max_y + min_y) / 2.0;
    double mid_xx = mov_x;
    double mid_yy = mov_y;
    if (mid_xx < EPS || mid_yy < EPS)
    {
    }
    else
    {
        mid_x = mid_xx;
        mid_y = mid_yy;
    }
    for (int i = 0; i < fileroot.size(); i++)//归一化
    {
        fileroot[i].x = fileroot[i].x - mid_x;
        fileroot[i].y = fileroot[i].y - mid_y;
    }
    geodown = min_z;
    day_max_z = max_z;
    path = path + bu;
//-------------------清空空面点云-------------------
//sort(fileroot.begin(), fileroot.end(), tsortup);
    int time_r = 0;
    while (time_r < fileroot.size())
    {
        int time_b = time_r;
        double limtime = fileroot[time_b].time + timegap * (double)86400000;
        while (fileroot[time_r].time < limtime)
        {
            time_r++;
            if (time_r >= fileroot.size())
                break;
        }
        string obj = ".obj";
        string time_name = to_string(num_time_name);
        num_day timemid;
        timemid.num_dd = num_time_name;
        timemid.day = fileroot[time_b].time;
        daylist.push_back(timemid);
        ofstream myout(path + time_name + obj);
        myout << "o " << num_time_name << "\n";
        num_time_name++;
        Polygo file2;//清空少点面
        vector<double> kongdown;//实面高度
        vector<floatz> full;
        double down = min_z;
//清空少点面
        clearempty(fileroot, kongdown, time_b, time_r, lim, file2, down, dz);
        if (file2.empty())
        {
            continue;
        }
        creatfull(file2, kongdown, full, dz);
        sort(full.begin(), full.end(), fsortup);
//---------------输出准备---------------------
        int read_num = 0;
        int num = 1;//顺序标记，输入点顺序
//------------切割z轴----------------
        down = geodown;//设定底层高度
        day_min_z = file2.front().z;
        day_max_z = daymaxz(file2, day_min_z, day_max_z, down);
        sort(file2.begin(), file2.end(), zsortup);
        read_num = 0;
        for (int full_read = 0; full_read < full.size(); full_read++)
        {
            Polygo zfloat;
            while (full[full_read].under - dz + EPS < file2[read_num].z && file2[read_num].z + EPS < full[full_read].up + dz)
            {
                zfloat.push_back(file2[read_num]);
                read_num++;
                if (read_num >= file2.size())
                    break;
            }
            if (zfloat.empty())
            {
                continue;
            }
            int z_num = zfloat.size();
            int beishu = (full[full_read].up - full[full_read].under) / dz;
// -----------查找边界---------------
            sort(zfloat.begin(), zfloat.end(), zsortup);
            Polygo Out = andrewScan(zfloat);
            for (int i = 0; i < Out.size(); i++)
            {
                Out[i].z = zfloat.front().z;
            }
            double ans = 0.0;
            Point mid;
            mid.x = Out.front().x;
            mid.y = Out.front().y;
            mid.z = Out.front().z;
            Out.push_back(mid);
            for (int i = 0; i < Out.size() - 1; i++) {
                ans += (Out[i].x * Out[i + 1].y - Out[i + 1].x * Out[i].y);
            }
            ans = ans / 2;
            double face = fabs(ans);
            if (z_num < (beishu + add_float) * min_point_num && face>max_face)
            {
                zfloat.clear();
                continue;
            }
            if (face < min_face)
            {
                zfloat.clear();
                continue;
            }
//------------按顺序输出点云---------------
            writeOBJ(Out, myout, num, dz, down);
        }
        myout.close();
    }
    cout << "建模结束" << endl;
//==============转gltf文件======================
    obj2gltf(path, num_time_name);
//=============写czml文件===================
//writeCZML(path, min_t, max_t, daylist, longitude, latitude, heightDeo);
//return 0;
}