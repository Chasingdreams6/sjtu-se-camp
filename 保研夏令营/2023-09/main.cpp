// @By chasingdreams
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
//#define DEBUG 1
#ifdef DEBUG
#define DPRINTF(format,args...) printf(format, ##args)
#else
#define DPRINTF(format,args...)
#endif
using namespace std;

/*
 * 要求图中最多的点数为N，
 * N可以设置为不超过内存限制的最大值；
 * 这里取N=1e3是因为large_data的范围是500
 * */
const int N = 1e3+5;
const int INF = 0x3f3f3f3f;

/*
 * 可在c++14下直接编译
 * 需要/data文件夹在程序当前目录
 * 示例编译方法："cmake.exe" --build .\ --target all -j 14
 * n,m,p,q 如题意
 * 采用邻接矩阵存图
 * weight[i][j]表示i->j的边权，-1表示没有边
 * dist[i][j]表示i->j的最短路径
 * indegree[i]表示i号点的入度
 * ps表示p个乘客
 * path记录randomWalk过程中的路径
 * edges[i]表示i号点的所有出边，事先排序加速randomWalk
 * Gen是随机数生成器
 * */
int n, m, q, p;
int weight[N][N];
int dist[N][N];
int indegree[N];
vector<int> ps;
vector<int> path;
vector<int> edges[N];
class MersenneTwister {
public:
    MersenneTwister() {
        index = 0;
        mt[0] = 2023; // seed
        for (int i = 1; i < 624; ++i) {
            mt[i] = (0x6C078965 * (mt[i-1] ^ (mt[i-1] >> 30)) + i)
                    & 0xFFFFFFFF;
        }
    }

    unsigned int generate() {
        if (index == 0) {
            for (int i = 0; i < 624; ++i) {
                unsigned int y =
                        ((mt[i] & 0x80000000) + (mt[(i + 1) % 624] & 0x7FFFFFFF))
                        & 0xFFFFFFFF;
                mt[i] = mt[(i + 397) % 624] ^ (y >> 1);
                if (y % 2 != 0) {
                    mt[i] ^= 0x9908B0DF;
                }
            }
        }

        unsigned int y = mt[index];
        y ^= (y >> 11);
        y ^= ((y << 7) & 0x9D2C5680);
        y ^= ((y << 15) & 0xEFC60000);
        y ^= (y >> 18);

        index = (index + 1) % 624;
        return y;
    }

    unsigned int mt[624];
    int index;
}Gen;


void init()
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            weight[i][j] = -1;
        }
    }
}

void BuildGraph()
{
    init();
    string fpath;
    cin >> fpath;
    ifstream fin;
    fin.open(fpath.c_str(), ios::in);
    fin >> n >> m >> q >> p;
    for (int i = 0; i < m; ++i) {
        int x, y, w;
        fin >> x >> y >> w;
        weight[x][y] = w; // add edge
        edges[x].emplace_back(y);
        indegree[y]++;
    }
    for (int i = 0; i < p; ++i) {
        int t;
        fin >> t;
        ps.emplace_back(t);
    }
    fin.close();
    int minD = indegree[0], minId = 0;
    for (int i = 1; i < n; ++i) {
        if (indegree[i] > minD) {
            minD = indegree[i];
            minId = i;
        }
    }
    cout << minId << " " << minD << endl;
}

void init_floyd()
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) dist[i][j] = 0;
            else {
                if (weight[i][j] != -1) dist[i][j] = weight[i][j];
                else dist[i][j] = INF;
            }
        }
    }
}

void debug()
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            DPRINTF("%d ", dist[i][j]);
        }
        DPRINTF("\n");
    }
}

/*
 * 使用floyd算法计算任意两点之间的最短路径
 * */
void floyd()
{
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
            }
        }
    }
    //debug();
}

void ComputePath()
{
    init_floyd();
    floyd();
    int ansx, ansy, answ = -1;
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            if (dist[x][y] > answ) {
                answ = dist[x][y];
                ansx = x;
                ansy = y;
            } else if (dist[x][y] == answ) {
                if (x < ansx) {
                    ansx = x;
                    ansy = y;
                } else if (x == ansx) {
                    if (y < ansy) {
                        ansy = y;
                    }
                }
            }
        }
    }
    cout << ansx << " " << ansy << " " << answ << endl;
}

/*
 * 将每个点的出边按照编号从小到大排序
 * */
void init_rw()
{
   for (int i = 0; i < n; ++i) {
       sort(edges[i].begin(), edges[i].end());
   }
}

/*
 * 模拟行走的过程
 * cur_done_pos表示当前已经送到目的地的人下标
 * cur_carry_pos表示已上车的人下标
 * cur_node表示当前在的点标号
 * */
void RandomWalk()
{
    init_rw();
    int cur_done_pos = -1, cur_carry_pos = 0, cur_node = 0;
    int cur_size = 0;
    long long totw = 0;
    while (cur_done_pos + 1 < ((int)ps.size())) { // 直到所有乘客被送到
        path.emplace_back(cur_node);
        // only one on each node
        if (cur_carry_pos < (int)ps.size() &&
            ps[cur_carry_pos] == cur_node &&
            cur_size < q) { // 尝试上车, 题目保证每个点最多只会上一个乘客
            cur_carry_pos++;
            cur_size++;
        }
        if (!cur_node) { // 在0号点下车
            cur_done_pos = cur_carry_pos - 1;
            cur_size = 0;
        }
        if (cur_done_pos == ps.size() - 1) { // 完成模拟
            break;
        }
        int nxt_node = edges[cur_node][Gen.generate() % (edges[cur_node].size())]; // cur_node->nxt_node
        totw += weight[cur_node][nxt_node];
        cur_node = nxt_node;
    }
    cout << totw << endl;
    ofstream fout;
    fout.open("path.txt", ios::out);
    for (auto cur : path) {
        fout << cur << endl;
    }
    fout.close();
}
/**
 * poi是用于排序的结构体，重载了<
 */
struct poi {
    int node;
    int cnt;
    bool operator < (const poi &t) const
    {
        if (cnt != t.cnt)
            return cnt > t.cnt;
        return node < t.node;
    }
}tmp[N];
int tmp_tot;
/*
 * path_cnt[x]表示x号点在路径中出现的次数
 * */
map<int, int> path_cnt;

void PathAnalyze()
{
    for (auto x : path) {
        if (!path_cnt.count(x)) {
            path_cnt[x] = 1;
        } else {
            path_cnt[x]++;
        }
    }
    for (auto pair : path_cnt) {
        tmp[++tmp_tot].cnt = pair.second;
        tmp[tmp_tot].node = pair.first;
    }
    sort(tmp + 1, tmp + tmp_tot + 1); // 对路径中的点进行排序
    for (int i = 1; i <= 3; ++i) {
        cout << tmp[i].node << " " << tmp[i].cnt << endl;
    }
}


int f[N][N];
/*
 * 初始化dp的变量
 *
 * */
void init_dp()
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            f[i][j] = INF;
        }
    }
    f[0][0] = 0;
}
/*
 * f[i][k]表示已经装了i个人(0<=i<=p)，使用了k的空间的最小代价
 * f[i][k]可以转移到 f[i+1][k+1]，表示直接从pos[i]到pos[i+1]，走最短路径
 * pos[i]直接到pos[i+1]并且途中经过0号点的情况被下面包含；
 * f[i][k]也可以转移到f[i+1][1]，表示先从pos[i]回0放人，再去pos[i+1];
 * 收敛条件：答案不会继续变小，停止dp
 * */
void loop() {
    bool flag = true;
    while (flag) {
        flag = false;
        for (int i = 0; i < p; ++i) { // 搭载pos=ps[i]这个位置的人，i=0时pos=0
            for (int k = 0; k <= q; ++k) { // 枚举当前已用的空间k
                if (f[i][k] != INF) {
                    int last, now; // last表示在运完第pos=ps[i-1]后车的位置
                    if (!i) {
                        last = 0;
                    } else {
                        last = ps[i - 1];
                    }
                    now = ps[i]; // now表示运第pos=ps[i]这个人的位置
                    if (k < q) { // last -> now
                        if (f[i][k] + dist[last][now] < f[i + 1][k + 1]) {
                            f[i + 1][k + 1] = f[i][k] + dist[last][now];
                            //DPRINTF("f[%d][%d] = %d, %d dir %d \n", i + 1, k + 1, f[i + 1][k + 1], last, now);
                            flag = true;
                        }
                    }
                    if (f[i + 1][1] > f[i][k] + dist[last][0] + dist[0][now]) { // last -> 0 -> now
                        f[i + 1][1] = f[i][k] + dist[last][0] + dist[0][now];
                        //DPRINTF("f[%d][%d] = %d, %d ret and go %d\n", i + 1, 1, f[i + 1][1], last, now);
                        flag = true;
                    }
                }
            }
        }
    }
    int ans = INF;
    // 统计答案，枚举最后所用的空间，并加上返程的距离
    for (int i = 1; i <= q; ++i) {
        //printf("usedspace = %d va=%d+%d=%d\n", i, f[p][i], dist[ps[q-1]][0], f[p][i] + dist[ps[q - 1]][0] );
        ans = min(ans, f[p][i] + dist[ps[p - 1]][0]);
    }
    cout << ans << endl;
}
void BestPath()
{
    init_dp();
    loop();
}

int main() {
    BuildGraph();
    ComputePath();
    RandomWalk();
    PathAnalyze();
    BestPath();
    return 0;
}
