#include "bits/stdc++.h"
using namespace std;
#define ll long long
#define u64 unsigned long long

int n, m, seed;

struct Edge {
    int dest;
    int weight;
};

struct Pivot{
    unordered_set<int> P;
    unordered_set<int> W;
};

struct Result {
    ll B;
    unordered_set<int> U;
};

class DataStructureD {
private:
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<pair<ll,int>>> heap;
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<pair<ll,int>>> heap_prepend;
    unordered_map<int, ll> best;
    int M;
    ll uppbound;

public:
    DataStructureD (int m, ll upp) : M(m), uppbound(upp) {}

    void insert(int v, ll val) {
        if (!best.count(v) or val < best[v]) {
            best[v] = val;
            heap.emplace(val,v);
        }
    }

    void Batch_Prepend(const vector<pair<int,ll>> &edges) {
        for (auto &edge : edges) {
            int v = edge.first;
            ll wt = edge.second;
            if (!best.count(v) or wt < best[v]) {
                best[v] = wt;
            }
        }

        if (edges.size() <= (size_t) M) for (auto &edge : edges) heap_prepend.emplace(edge.second, edge.first);
        else {
            for (int i = 0; i < edges.size(); i += (M/2)) {
                int curr = min(M/2, (int) edges.size()-i);
                for (int j = 0; j < curr; j++) {
                    heap_prepend.emplace(edges[i+j].second, edges[i+j].first);
                }
            }
        }
    }

    bool isEmpty() {
        return heap.empty() and heap_prepend.empty();
    }

    Result pull() {
        if (heap.empty() and heap_prepend.empty()) return {uppbound, {}};
        
        unordered_set<int> Si;
        ll Bi = uppbound;
        if (!heap_prepend.empty()) Bi = heap_prepend.top().first;
        while (!heap_prepend.empty() and Si.size() < (size_t) M) {
            int v = heap_prepend.top().second;
            ll wt = heap_prepend.top().first;
            heap_prepend.pop();
            
            if (best.count(v) and best[v] == wt) {
                Si.insert(v);
                best.erase(v);
            }
        }

        if (Si.size() >= (size_t) M) return {Bi,Si};

        if (!heap.empty()) Bi = min(Bi, heap.top().first);
        while (!heap.empty() and Si.size() < (size_t) M) {
            int v = heap.top().second;
            ll wt = heap.top().first;
            heap.pop();
            if (best.count(v) and best[v] == wt) {
                Si.insert(v);
                best.erase(v);
            }
        }


        return {Bi,Si};
    }
};

Pivot FindPivot(int k, const vector<vector<Edge>> &edges, ll B, const unordered_set<int> &S, vector<ll> &dist, int p_limit) {
    // cout << "S: ";
    // for (int s : S) cout << s << " ";
    // cout << "\n";
    vector<int> cand;
    for (int u : S)
        if (dist[u] < B) cand.push_back(u);

    unordered_set<int> P;
    if (cand.empty()) {
        int cnt = min((int)S.size(), max(1,p_limit));
        auto it = S.begin();
        while (cnt-- && it != S.end()) P.insert(*it++);
    } else {
        sort(cand.begin(), cand.end(),
             [&](int a, int b){ return dist[a] < dist[b]; });
        for (int i = 0; i < min((int)cand.size(), p_limit); i++)
            P.insert(cand[i]);
    }

    unordered_set<int> W = P.empty() ? S : P;
    unordered_set<int> frontier = W;

    for (int step = 0; step < max(1,k); step++) {
        unordered_set<int> nxt;
        for (int u : frontier) {
            if (dist[u] >= B) continue;
            for (auto &e : edges[u]) {
                ll nd = dist[u] + e.weight;
                if (nd < B && !W.count(e.dest)) {
                    W.insert(e.dest);
                    nxt.insert(e.dest);
                }
            }
        }
        if (nxt.empty()) break;
        frontier.swap(nxt);
    }

    if (P.empty() && !S.empty()) {P.insert(*S.begin()); return {P,W};}

    // vector<int> F(n, -1);
    // vector<int> subtree_size(n,0);
    // for (int u : W) {
    //     for (const Edge& edge : edges[u]) {
    //         int v = edge.dest;
    //         int wt = edge.weight;

    //         if (dist[v] == dist[u] + wt) {F[v] = u; subtree_size[F[v]]++;}
    //     }
    // }

    // for (int u : S) {
    //     if (subtree_size[u] >= k) P.insert(u);
    // }

    return {P,W};
}

Result BaseCase(int k, const vector<vector<Edge>> &edges, ll B, const unordered_set<int> &S, vector<ll> &dist) {
    int x = *S.begin();
    unordered_set<int> U0;
    U0.insert(x);
    priority_queue<pair<ll,int>, vector<pair<ll, int>>, greater<pair<ll,int>>> H;
    H.emplace(dist[x],x);

    while (!H.empty() and U0.size() < (size_t) k+1) {
        pair<ll,int> u = H.top();
        H.pop();
        if (u.first > dist[u.second]) continue;
        U0.insert(u.second);

        for (const Edge &edge : edges[u.second]) {
            int v = edge.dest;
            int wt = edge.weight;

            if (dist[u.second] + wt <= dist[v] and dist[u.second] + wt < B) {
                dist[v] = dist[u.second] + wt;
                H.emplace(dist[v],v);
            }
        }
    }
    // cout << "BASE CASE B = " << B << ": ";
    // for (int u : U0) cout << u << " ";
    // cout << "\n";
    if (U0.size() <= (size_t) k) return {B, U0};
    else {
        ll B_max = -1;
        for (int v : U0) {
            if (dist[v] > B_max) B_max = dist[v];
        }

        unordered_set<int> U;
        for (int v : U0) {
            if (dist[v] < B_max) U.insert(v);
        }
        return {B_max,U};
    }
}

Result BMSSP(int l, const vector<vector<Edge>> &edges, ll B, const unordered_set<int> &S, vector<ll> &dist) {
    // cout << "l: " << l << "\n";
    int k = max(1, (int) pow(log2(n), 1.0/3.0));
    int t = max(1, (int) pow(log2(n), 2.0/3.0));
    int M = (int) pow(2, max(0,(l-1)*t));
    int p_limit = max(1, (int) pow(2, min(10,t)));
    if (l == 0) {
        return BaseCase(k,edges,B,S,dist);
    }

    Pivot pivot = FindPivot(k,edges,B,S,dist, p_limit);
    // cout << "P: ";
    // for (int p : pivot.P) cout << p << " ";
    // cout << "\nW: ";
    // for (int w : pivot.W) cout << w << " ";
    // cout << "\n";

    DataStructureD D(M,B);
    
    for (int x : pivot.P) D.insert(x, dist[x]);

    ll B0 = B;
    if (!pivot.P.empty() and *max_element(pivot.P.begin(), pivot.P.end()) != 0) {
        for (int x : pivot.P) B0 = min(B0, dist[x]);
    }
    
    unordered_set<int> U;
    while (U.size() < (size_t) k*(pow(2,l*t)) and !D.isEmpty()) {
        Result ans = D.pull();
        // cout << "B: " << ans.B << "\n";
        // cout << "U: ";
        // for (int u : ans.U) cout << u << " ";
        // cout << "\n";
        Result ans1 = BMSSP(l-1,edges, ans.B, ans.U, dist);
        B0 = min(B0,ans1.B);
        
        vector<pair<int,ll>> K;
        for (int u : ans1.U) {
            U.insert(u);
            for (const Edge &edge : edges[u]) {
                int v = edge.dest;
                int wt = edge.weight;
                if (dist[u] != LLONG_MAX and dist[u] + wt <= dist[v]) {
                    dist[v] = dist[u] + wt;
                    if (dist[u] + wt < B and dist[u] + wt >= ans.B) D.insert(v, dist[u] + wt);
                    else if (dist[u] + wt < ans.B and dist[u] + wt >= B0) K.push_back({v,dist[u] + wt});
                }
            }
        }
        for (int x : ans.U) {
            if (dist[x] < ans.B and dist[x] >= B0) K.push_back({x, dist[x]});
        }
        if (!K.empty()) D.Batch_Prepend(K);
    }

    B0 = min(B, B0);
    for (int x : pivot.W) {
        if (dist[x] < B0) U.insert(x); 
    }
    // cout << "Dist: ";
    // for (ll d : dist) cout << d << " ";
    // cout << "\n";
    return {B0,U};
}

void dijkstra(int src, const vector<vector<Edge>>& adj, vector<ll>& dist) {
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<>> pq;
    pq.emplace(0, src);
    while (!pq.empty()) {
        auto [d,u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        for (auto &e : adj[u]) {
            if (d + e.weight < dist[e.dest]) {
                dist[e.dest] = d + e.weight;
                pq.emplace(dist[e.dest], e.dest);
            }
        }
    }
}

int main() {
    cin >> n >> m >> seed;
    int src = 0;

    vector<vector<Edge>> adj(n);

    mt19937 rng(seed);
    uniform_int_distribution<int> W(1,100);

    // Backbone
    for (int i = 1; i < n; i++) {
        uniform_int_distribution<int> U(0,i-1);
        adj[U(rng)].push_back({i, W(rng)});
    }

    // // Remaining edges
    // for (int i = 0; i < m-(n-1); i++) {
    //     int u = rng()%n, v = rng()%n;
    //     if (u != v) adj[u].push_back({v, W(rng)});
    // }

    vector<ll> d_dij(n, LLONG_MAX);
    d_dij[src] = 0;

    auto t1 = chrono::high_resolution_clock::now();
    dijkstra(src, adj, d_dij);
    auto t2 = chrono::high_resolution_clock::now();

    double time_dij = chrono::duration<double>(t2 - t1).count();

    vector<ll> d_bm(n, LLONG_MAX);
    d_bm[src] = 0;

    int t = max(1, (int)round(pow(log(max(3,n)), 2.0/3.0)));
    int l = max(1, (int)round(log(max(3,n)) / t));

    auto t3 = chrono::high_resolution_clock::now();
    BMSSP(l, adj, LLONG_MAX, {src}, d_bm);
    auto t4 = chrono::high_resolution_clock::now();

    double time_bm = chrono::duration<double>(t4 - t3).count();

    ll max_diff = 0;
    int common_reachable = 0;

    for (int i = 0; i < n; i++) {
        if (d_dij[i] < LLONG_MAX && d_bm[i] < LLONG_MAX) {
            max_diff = max(max_diff, llabs(d_dij[i] - d_bm[i]));
            common_reachable++;
        }
    }

    cout << "Dijkstra time: " << time_dij << " s\n";
    cout << "BMSSP time   : " << time_bm << " s\n";
    cout << "Max diff     : " << max_diff << "\n";
    cout << "Common reachable: " << common_reachable << "\n";

    return 0;
}
