#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
#define int long long
#define vi vector<int>
#define pb push_back
#define mp make_pair
#define ff first
#define ss second
#define inar(arr,n) for(int i=0; i<n; i++) cin>>arr[i];
#define suma(arr,n) for(int i=0; i<n; i++) sum+=arr[i];
#define for0(i,n) for(int i=0; i<n; i++)
#define for1(i,n) for(int i=1; i<=n; i++)
#define forr(i,n) for(int i=n-1; i>=0; i--)
#define rof(i,a,b) for(int i=a;i<=b;i++)
#define all(v) v.begin(), v.end()
#define rall(v) v.rbegin(), v.rend()
#define print(v) cout<<v<<"\n"
#define fast ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0)
const int mod = 1e9 + 7;
ll gcd (ll a, ll b) {return ( a ? gcd(b % a, a) : b );}
ll power(ll a, ll n) {ll p = 1; while (n > 0) {if (n % 2) {p = p * a;} n >>= 1; a *= a;} return p;}
ll power(ll a, ll n, ll mod) {ll p = 1; while (n > 0) {if (n % 2) {p = p * a; p %= mod;} n >>= 1; a *= a; a %= mod;} return p % mod;}

signed main() {
  fast;
#ifndef ONLINE_JUDGE
  freopen("input.txt", "r", stdin);
  freopen("output.txt", "w", stdout);
#endif

  // int k = 18;
  // for (int x = 0; x <= k; x++) {
  //   if ((power(x, 2) + 2 * x + 3) % k == 0) {
  //     cout << x << " ";
  //   }
  // }

  vector<int> a1, a2;
  a1 = {1};
  a2 = {3, 4};
  for (int i : a1) {
    for (int j : a2) {
      cout << ((8 * j * -1 + 9 * i) + 36) % 18 << " ";
    }
  }

  // cout << power(331, 2, 495);
  return 0;
}