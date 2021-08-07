#ifndef NAGISS_LIBRARY_HPP
#define NAGISS_LIBRARY_HPP
#include<iostream>
#include<iomanip>
#include<vector>
#include<set>
#include<map>
#include<unordered_set>
#include<unordered_map>
#include<algorithm>
#include<numeric>
#include<limits>
#include<bitset>
#include<functional>
#include<type_traits>
#include<queue>
#include<stack>
#include<array>
#include<random>
#include<utility>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<string>
#include<sstream>
#include<chrono>
#include<climits>
#ifdef _MSC_VER
#include<intrin0.h>
#endif

#ifdef __GNUC__
//#pragma GCC target("avx2")
//#pragma GCC target("sse4")
//#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
//#pragma GCC optimize("O3")
//#pragma GCC optimize("Ofast")
//#pragma GCC optimize("unroll-loops")
#endif

// ========================== macroes ==========================

#define rep(i,n) for(ll (i)=0; (i)<(n); (i)++)
#define rep1(i,n) for(ll (i)=1; (i)<=(n); (i)++)
#define rep3(i,s,n) for(ll (i)=(s); (i)<(n); (i)++)

#define NDEBUG

#ifndef NDEBUG
#define ASSERT(expr, ...) \
		do { \
			if(!(expr)){ \
				printf("%s(%d): Assertion failed.\n", __FILE__, __LINE__); \
				printf(__VA_ARGS__); \
				abort(); \
			} \
		} while (false)
#else
#define ASSERT(...)
#endif

#define ASSERT_RANGE(value, left, right) \
	ASSERT((left <= value) && (value < right), \
		"`%s` (%d) is out of range [%d, %d)", #value, value, left, right)

#define CHECK(var) do{ std::cout << #var << '=' << var << endl; } while (false)

// ========================== utils ==========================

using namespace std;
using ll = long long;
constexpr double PI = 3.1415926535897932;

template<class T, class S> inline bool chmin(T& m, S q) {
	if (m > q) { m = q; return true; }
	else return false;
}

template<class T, class S> inline bool chmax(T& m, const S q) {
	if (m < q) { m = q; return true; }
	else return false;
}

// クリッピング  // clamp (C++17) と等価
template<class T> inline T clipped(const T& v, const T& low, const T& high) {
	return min(max(v, low), high);
}

// 2 次元ベクトル
template<typename T> struct Vec2 {
	/*
	y 軸正は下方向
	x 軸正は右方向
	回転は時計回りが正（y 軸正を上と考えると反時計回りになる）
	*/
	T y, x;
	constexpr inline Vec2() = default;
	constexpr Vec2(const T& arg_y, const T& arg_x) : y(arg_y), x(arg_x) {}
	inline Vec2(const Vec2&) = default;  // コピー
	inline Vec2(Vec2&&) = default;  // ムーブ
	inline Vec2& operator=(const Vec2&) = default;  // 代入
	inline Vec2& operator=(Vec2&&) = default;  // ムーブ代入
	template<typename S> constexpr inline Vec2(const Vec2<S>& v) : y((T)v.y), x((T)v.x) {}
	inline Vec2 operator+(const Vec2& rhs) const {
		return Vec2(y + rhs.y, x + rhs.x);
	}
	inline Vec2 operator+(const T& rhs) const {
		return Vec2(y + rhs, x + rhs);
	}
	inline Vec2 operator-(const Vec2& rhs) const {
		return Vec2(y - rhs.y, x - rhs.x);
	}
	template<typename S> inline Vec2 operator*(const S& rhs) const {
		return Vec2(y * rhs, x * rhs);
	}
	inline Vec2 operator*(const Vec2& rhs) const {  // x + yj とみなす
		return Vec2(x * rhs.y + y * rhs.x, x * rhs.x - y * rhs.y);
	}
	template<typename S> inline Vec2 operator/(const S& rhs) const {
		ASSERT(rhs != 0.0, "Zero division!");
		return Vec2(y / rhs, x / rhs);
	}
	inline Vec2 operator/(const Vec2& rhs) const {  // x + yj とみなす
		return (*this) * rhs.inv();
	}
	inline Vec2& operator+=(const Vec2& rhs) {
		y += rhs.y;
		x += rhs.x;
		return *this;
	}
	inline Vec2& operator-=(const Vec2& rhs) {
		y -= rhs.y;
		x -= rhs.x;
		return *this;
	}
	template<typename S> inline Vec2& operator*=(const S& rhs) const {
		y *= rhs;
		x *= rhs;
		return *this;
	}
	inline Vec2& operator*=(const Vec2& rhs) {
		*this = (*this) * rhs;
		return *this;
	}
	inline Vec2& operator/=(const Vec2& rhs) {
		*this = (*this) / rhs;
		return *this;
	}
	inline bool operator!=(const Vec2& rhs) const {
		return x != rhs.x || y != rhs.y;
	}
	inline bool operator==(const Vec2& rhs) const {
		return x == rhs.x && y == rhs.y;
	}
	inline void rotate(const double& rad) {
		*this = rotated(rad);
	}
	inline Vec2<double> rotated(const double& rad) const {
		return (*this) * rotation(rad);
	}
	static inline Vec2<double> rotation(const double& rad) {
		return Vec2(sin(rad), cos(rad));
	}
	static inline Vec2<double> rotation_deg(const double& deg) {
		return rotation(PI * deg / 180.0);
	}
	inline Vec2<double> rounded() const {
		return Vec2<double>(round(y), round(x));
	}
	inline Vec2<double> inv() const {  // x + yj とみなす
		const double norm_sq = l2_norm_square();
		ASSERT(norm_sq != 0.0, "Zero division!");
		return Vec2(-y / norm_sq, x / norm_sq);
	}
	inline double l2_norm() const {
		return sqrt(x * x + y * y);
	}
	inline double l2_norm_square() const {
		return x * x + y * y;
	}
	inline T l1_norm() const {
		return std::abs(x) + std::abs(y);
	}
	inline double abs() const {
		return l2_norm();
	}
	inline double phase() const {  // [-PI, PI) のはず
		return atan2(y, x);
	}
	inline double phase_deg() const {  // [-180, 180) のはず
		return phase() / PI * 180.0;
	}
};
template<typename T, typename S> inline Vec2<T> operator*(const S& lhs, const Vec2<T>& rhs) {
	return rhs * lhs;
}
template<typename T> ostream& operator<<(ostream& os, const Vec2<T>& vec) {
	os << vec.y << ' ' << vec.x;
	return os;
}

// 乱数
struct Random {
	using ull = unsigned long long;
	ull seed;
	inline Random(ull aSeed) : seed(aSeed) {
		ASSERT(seed != 0ull, "Seed should not be 0.");
	}
	const inline ull& next() {
		seed ^= seed << 9;
		seed ^= seed >> 7;
		return seed;
	}
	// (0.0, 1.0)
	inline double random() {
		return (double)next() / (double)ULLONG_MAX;
	}
	// [0, right)
	inline int randint(const int right) {
		return next() % (ull)right;
	}
	// [left, right)
	inline int randint(const int left, const int right) {
		return next() % (ull)(right - left) + left;
	}
};


// キュー
template<class T, int max_size> struct Queue {
	array<T, max_size> data;
	int left, right;
	inline Queue() : data(), left(0), right(0) {}
	inline Queue(initializer_list<T> init) :
		data(init.begin(), init.end()), left(0), right(init.size()) {}

	inline bool empty() const {
		return left == right;
	}
	inline void push(const T& value) {
		data[right] = value;
		right++;
	}
	inline void pop() {
		left++;
	}
	const inline T& front() const {
		return data[left];
	}
	template <class... Args> inline void emplace(const Args&... args) {
		data[right] = T(args...);
		right++;
	}
	inline void clear() {
		left = 0;
		right = 0;
	}
	inline int size() const {
		return right - left;
	}
};


// スタック  // コンストラクタ呼ぶタイミングとかが考えられてなくて良くない
template<class T, int max_size> struct Stack {
	array<T, max_size> data;
	int right;
	inline Stack() : data(), right(0) {}
	inline Stack(const int n) : data(), right(0) { resize(n); }
	inline Stack(const int n, const T& val) : data(), right(0) { resize(n, val); }
	inline Stack(initializer_list<T> init) :
		data(init.begin(), init.end()), right(init.size()) {}  // これ長さが最大じゃないと動かない
	inline Stack(const Stack& rhs) : data(), right(rhs.right) {  // コピー
		for (int i = 0; i < right; i++) {
			data[i] = rhs.data[i];
		}
	}
	Stack& operator=(const Stack& rhs) {
		right = rhs.right;
		for (int i = 0; i < right; i++) {
			data[i] = rhs.data[i];
		}
		return *this;
	}
	Stack& operator=(const vector<T>& rhs) {
		right = (int)rhs.size();
		ASSERT(right <= max_size, "too big vector");
		for (int i = 0; i < right; i++) {
			data[i] = rhs[i];
		}
		return *this;
	}
	Stack& operator=(Stack&&) = default;
	inline bool empty() const {
		return 0 == right;
	}
	inline void push(const T& value) {
		ASSERT_RANGE(right, 0, max_size);
		data[right] = value;
		right++;
	}
	inline T pop() {
		right--;
		ASSERT_RANGE(right, 0, max_size);
		return data[right];
	}
	const inline T& top() const {
		return data[right - 1];
	}
	template <class... Args> inline void emplace(const Args&... args) {
		ASSERT_RANGE(right, 0, max_size);
		data[right] = T(args...);
		right++;
	}
	inline void clear() {
		right = 0;
	}
	inline void insert(const int& idx, const T& value) {
		ASSERT_RANGE(idx, 0, right + 1);
		ASSERT_RANGE(right, 0, max_size);
		int i = right;
		right++;
		while (i != idx) {
			data[i] = data[i - 1];
			i--;
		}
		data[idx] = value;
	}
	inline void del(const int& idx) {
		ASSERT_RANGE(idx, 0, right);
		right--;
		for (int i = idx; i < right; i++) {
			data[i] = data[i + 1];
		}
	}
	inline int index(const T& value) const {
		for (int i = 0; i < right; i++) {
			if (value == data[i]) return i;
		}
		return -1;
	}
	inline void remove(const T& value) {
		int idx = index(value);
		ASSERT(idx != -1, "not contain the value.");
		del(idx);
	}
	inline void resize(const int& sz) {
		ASSERT_RANGE(sz, 0, max_size + 1);
		for (; right < sz; right++) {
			data[right].~T();
			new(&data[right]) T();
		}
		right = sz;
	}
	inline void resize(const int& sz, const T& fill_value) {
		ASSERT_RANGE(sz, 0, max_size + 1);
		for (; right < sz; right++) {
			data[right].~T();
			new(&data[right]) T(fill_value);
		}
		right = sz;
	}
	inline int size() const {
		return right;
	}
	inline T& operator[](const int n) {
		ASSERT_RANGE(n, 0, right);
		return data[n];
	}
	inline const T& operator[](const int n) const {
		ASSERT_RANGE(n, 0, right);
		return data[n];
	}
	inline T* begin() {
		return (T*)data.data();
	}
	inline const T* begin() const {
		return (const T*)data.data();
	}
	inline T* end() {
		return (T*)data.data() + right;
	}
	inline const T* end() const {
		return (const T*)data.data() + right;
	}
	inline T& front() {
		ASSERT(right > 0, "no data.");
		return data[0];
	}
	const inline T& front() const {
		ASSERT(right > 0, "no data.");
		return data[0];
	}
	inline T& back() {
		ASSERT(right > 0, "no data.");
		return data[right - 1];
	}
	const inline T& back() const {
		ASSERT(right > 0, "no data.");
		return data[right - 1];
	}
	inline bool contains(const T& value) const {
		for (const auto& dat : *this) {
			if (value == dat) return true;
		}
		return false;
	}
	inline vector<T> ToVector() {
		return vector<T>(begin(), end());
	}
};


// 時間 (秒)
inline double time() {
	return static_cast<double>(chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now().time_since_epoch()).count()) * 1e-9;
}


// 重複除去
template<typename T> inline void deduplicate(vector<T>& vec) {
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
}


template<typename T> inline int search_sorted(const vector<T>& vec, const T& a) {
	return lower_bound(vec.begin(), vec.end(), a) - vec.begin();
}


// popcount  // SSE 4.2 を使うべき
inline int popcount(const unsigned int& x) {
#ifdef _MSC_VER
	return (int)__popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}
inline int popcount(const unsigned long long& x) {
#ifdef _MSC_VER
	return (int)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
}

// x >> n & 1 が 1 になる最小の n ( x==0 は未定義 )
inline int CountRightZero(const unsigned int& x) {
#ifdef _MSC_VER
	unsigned long r;
	_BitScanForward(&r, x);
	return (int)r;
#else
	return __builtin_ctz(x);
#endif
}
inline int CountRightZero(const unsigned long long& x) {
#ifdef _MSC_VER
	unsigned long r;
	_BitScanForward64(&r, x);
	return (int)r;
#else
	return __builtin_ctzll(x);
#endif
}

#endif  // NAGISS_LIBRARY_HPP

#ifdef _MSC_VER
inline unsigned int __builtin_clz(const unsigned int& x) { unsigned long r; _BitScanReverse(&r, x); return 31 - r; }
inline unsigned long long __builtin_clzll(const unsigned long long& x) { unsigned long r; _BitScanReverse64(&r, x); return 63 - r; }
#endif
#include<cassert>
#pragma warning( disable : 4146 )
/*
iwi 先生の radix heap (https://github.com/iwiwi/radix-heap)
The MIT License (MIT)
Copyright (c) 2015 Takuya Akiba
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
namespace radix_heap {
namespace internal {
template<bool Is64bit> class find_bucket_impl;

template<>
class find_bucket_impl<false> {
public:
	static inline constexpr size_t find_bucket(uint32_t x, uint32_t last) {
		return x == last ? 0 : 32 - __builtin_clz(x ^ last);
	}
};

template<>
class find_bucket_impl<true> {
public:
	static inline constexpr size_t find_bucket(uint64_t x, uint64_t last) {
		return x == last ? 0 : 64 - __builtin_clzll(x ^ last);
	}
};

template<typename T>
inline constexpr size_t find_bucket(T x, T last) {
	return find_bucket_impl<sizeof(T) == 8>::find_bucket(x, last);
}

template<typename KeyType, bool IsSigned> class encoder_impl_integer;

template<typename KeyType>
class encoder_impl_integer<KeyType, false> {
public:
	typedef KeyType key_type;
	typedef KeyType unsigned_key_type;

	inline static constexpr unsigned_key_type encode(key_type x) {
		return x;
	}

	inline static constexpr key_type decode(unsigned_key_type x) {
		return x;
	}
};

template<typename KeyType>
class encoder_impl_integer<KeyType, true> {
public:
	typedef KeyType key_type;
	typedef typename std::make_unsigned<KeyType>::type unsigned_key_type;

	inline static constexpr unsigned_key_type encode(key_type x) {
		return static_cast<unsigned_key_type>(x) ^
			(unsigned_key_type(1) << unsigned_key_type(std::numeric_limits<unsigned_key_type>::digits - 1));
	}

	inline static constexpr key_type decode(unsigned_key_type x) {
		return static_cast<key_type>
			(x ^ (unsigned_key_type(1) << (std::numeric_limits<unsigned_key_type>::digits - 1)));
	}
};

template<typename KeyType, typename UnsignedKeyType>
class encoder_impl_decimal {
public:
	typedef KeyType key_type;
	typedef UnsignedKeyType unsigned_key_type;

	inline static constexpr unsigned_key_type encode(key_type x) {
		return raw_cast<key_type, unsigned_key_type>(x) ^
			((-(raw_cast<key_type, unsigned_key_type>(x) >> (std::numeric_limits<unsigned_key_type>::digits - 1))) |
				(unsigned_key_type(1) << (std::numeric_limits<unsigned_key_type>::digits - 1)));
	}

	inline static constexpr key_type decode(unsigned_key_type x) {
		return raw_cast<unsigned_key_type, key_type>
			(x ^ (((x >> (std::numeric_limits<unsigned_key_type>::digits - 1)) - 1) |
				(unsigned_key_type(1) << (std::numeric_limits<unsigned_key_type>::digits - 1))));
	}

private:
	template<typename T, typename U>
	union raw_cast {
	public:
		constexpr raw_cast(T t) : t_(t) {}
		operator U() const { return u_; }

	private:
		T t_;
		U u_;
	};
};

template<typename KeyType>
class encoder : public encoder_impl_integer<KeyType, std::is_signed<KeyType>::value> {};
template<>
class encoder<float> : public encoder_impl_decimal<float, uint32_t> {};
template<>
class encoder<double> : public encoder_impl_decimal<double, uint64_t> {};
}  // namespace internal

template<typename KeyType, typename EncoderType = internal::encoder<KeyType>>
class radix_heap {
public:
	typedef KeyType key_type;
	typedef EncoderType encoder_type;
	typedef typename encoder_type::unsigned_key_type unsigned_key_type;

	radix_heap() : size_(0), last_(), buckets_() {
		buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
	}

	void push(key_type key) {
		const unsigned_key_type x = encoder_type::encode(key);
		assert(last_ <= x);
		++size_;
		const size_t k = internal::find_bucket(x, last_);
		buckets_[k].emplace_back(x);
		buckets_min_[k] = std::min(buckets_min_[k], x);
	}

	key_type top() {
		pull();
		return encoder_type::decode(last_);
	}

	void pop() {
		pull();
		buckets_[0].pop_back();
		--size_;
	}

	size_t size() const {
		return size_;
	}

	bool empty() const {
		return size_ == 0;
	}

	void clear() {
		size_ = 0;
		last_ = key_type();
		for (auto& b : buckets_) b.clear();
		buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
	}

	void swap(radix_heap<KeyType, EncoderType>& a) {
		std::swap(size_, a.size_);
		std::swap(last_, a.last_);
		buckets_.swap(a.buckets_);
		buckets_min_.swap(a.buckets_min_);
	}

private:
	size_t size_;
	unsigned_key_type last_;
	std::array<std::vector<unsigned_key_type>,
		std::numeric_limits<unsigned_key_type>::digits + 1> buckets_;
	std::array<unsigned_key_type,
		std::numeric_limits<unsigned_key_type>::digits + 1> buckets_min_;

	void pull() {
		assert(size_ > 0);
		if (!buckets_[0].empty()) return;

		size_t i;
		for (i = 1; buckets_[i].empty(); ++i);
		last_ = buckets_min_[i];

		for (unsigned_key_type x : buckets_[i]) {
			const size_t k = internal::find_bucket(x, last_);
			buckets_[k].emplace_back(x);
			buckets_min_[k] = std::min(buckets_min_[k], x);
		}
		buckets_[i].clear();
		buckets_min_[i] = std::numeric_limits<unsigned_key_type>::max();
	}
};

template<typename KeyType, typename ValueType, typename EncoderType = internal::encoder<KeyType>>
class pair_radix_heap {
public:
	typedef KeyType key_type;
	typedef ValueType value_type;
	typedef EncoderType encoder_type;
	typedef typename encoder_type::unsigned_key_type unsigned_key_type;

	pair_radix_heap() : size_(0), last_(), buckets_() {
		buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
	}

	void push(key_type key, const value_type& value) {
		const unsigned_key_type x = encoder_type::encode(key);
		assert(last_ <= x);
		++size_;
		const size_t k = internal::find_bucket(x, last_);
		buckets_[k].emplace_back(x, value);
		buckets_min_[k] = std::min(buckets_min_[k], x);
	}

	void push(key_type key, value_type&& value) {
		const unsigned_key_type x = encoder_type::encode(key);
		assert(last_ <= x);
		++size_;
		const size_t k = internal::find_bucket(x, last_);
		buckets_[k].emplace_back(x, std::move(value));
		buckets_min_[k] = std::min(buckets_min_[k], x);
	}

	template <class... Args>
	void emplace(key_type key, Args&&... args) {
		const unsigned_key_type x = encoder_type::encode(key);
		assert(last_ <= x);
		++size_;
		const size_t k = internal::find_bucket(x, last_);
		buckets_[k].emplace_back(std::piecewise_construct,
			std::forward_as_tuple(x), std::forward_as_tuple(args...));
		buckets_min_[k] = std::min(buckets_min_[k], x);
	}

	key_type top_key() {
		pull();
		return encoder_type::decode(last_);
	}

	value_type& top_value() {
		pull();
		return buckets_[0].back().second;
	}

	void pop() {
		pull();
		buckets_[0].pop_back();
		--size_;
	}

	size_t size() const {
		return size_;
	}

	bool empty() const {
		return size_ == 0;
	}

	void clear() {
		size_ = 0;
		last_ = key_type();
		for (auto& b : buckets_) b.clear();
		buckets_min_.fill(std::numeric_limits<unsigned_key_type>::max());
	}

	void swap(pair_radix_heap<KeyType, ValueType, EncoderType>& a) {
		std::swap(size_, a.size_);
		std::swap(last_, a.last_);
		buckets_.swap(a.buckets_);
		buckets_min_.swap(a.buckets_min_);
	}

private:
	size_t size_;
	unsigned_key_type last_;
	std::array<std::vector<std::pair<unsigned_key_type, value_type>>,
		std::numeric_limits<unsigned_key_type>::digits + 1> buckets_;
	std::array<unsigned_key_type,
		std::numeric_limits<unsigned_key_type>::digits + 1> buckets_min_;

	void pull() {
		assert(size_ > 0);
		if (!buckets_[0].empty()) return;

		size_t i;
		for (i = 1; buckets_[i].empty(); ++i);
		last_ = buckets_min_[i];

		for (size_t j = 0; j < buckets_[i].size(); ++j) {
			const unsigned_key_type x = buckets_[i][j].first;
			const size_t k = internal::find_bucket(x, last_);
			buckets_[k].emplace_back(std::move(buckets_[i][j]));
			buckets_min_[k] = std::min(buckets_min_[k], x);
		}
		buckets_[i].clear();
		buckets_min_[i] = std::numeric_limits<unsigned_key_type>::max();
	}
};
}  // namespace radix_heap

auto rng = Random(42);

using Point = Vec2<signed char>;
int N, sy, sx;
double t0;
array<array<signed char, 71>, 71> board;

constexpr int MAX_N_CROSSINGS = 500;  // 足りるか？
auto n_roads = 0;
auto crossings = Stack<Point, MAX_N_CROSSINGS>();
auto visibles = array<Stack<Point, 69 * 2 - 1>, MAX_N_CROSSINGS>();
auto distances = array<array<short, MAX_N_CROSSINGS>, MAX_N_CROSSINGS>();
auto point_to_crossing_idx = array<array<short, 71>, 71>();

namespace dijkstra {

constexpr auto inf = numeric_limits<short>::max();
auto dists = array<array<short, 71>, 71>();
auto from = array<array<Point, 71>, 71>();
auto crossing_dists = array<short, MAX_N_CROSSINGS>();

void Dijkstra(const short& start_idx_crossing) {
	for (auto&& d : dists) for (auto&& dd : d) dd = inf;
	const auto& start = crossings[start_idx_crossing];
	from[start.y][start.x] = Point(-1, -1);
	dists[start.y][start.x] = 0;
	auto q = radix_heap::pair_radix_heap<short, Point>();
	q.push(0, start);
	while (q.size()) {
		const auto dist_v = q.top_key();
		const auto v = q.top_value();
		q.pop();
		if (dist_v != dists[v.y][v.x]) continue;
		for (const auto& dyx : { Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1) }) {
			const auto u = v + dyx;
			const auto& cost = board[u.y][u.x];
			if (cost == -1) continue;
			const auto dist_u = dist_v + cost;
			if (dist_u < dists[u.y][u.x]) {
				dists[u.y][u.x] = dist_u;
				from[u.y][u.x] = v;
				q.push(dist_u, u);
			}
		}
	}

	rep(idx_crossings, crossings.size()) {
		const auto& crossing = crossings[idx_crossings];
		crossing_dists[idx_crossings] = dists[crossing.y][crossing.x];
	}
}
}  // namespace dijkstra

using dijkstra::Dijkstra;

void Solve() {

	// input
	{
		for (auto&& b : board) for (auto&& bb : b) bb = -1;
		cin >> N >> sy >> sx;
		sy++; sx++;
		rep(y, N) {
			string s;
			cin >> s;
			rep(x, N) {
				char c = s[x];
				if (c != '#') board[y + 1][x + 1] = c - '0';
			}
		}
	}


	t0 = time();
	// 前処理
	// 交差点を列挙
	{
		crossings.emplace(sy, sx);
		rep1(y, N) rep1(x, N) {
			point_to_crossing_idx[y][x] = -1;
			auto p = Point(y, x);
			if (board[y][x] == -1) {
				continue;
			}
			n_roads++;
			if (y == sy && x == sx) {
				point_to_crossing_idx[y][x] = 0;
				continue;
			}
			int road_cnt = 0;
			for (const auto& dyx : { Point(0, 1), Point(0, -1) }) {
				auto uyx = p + dyx;
				if (board[uyx.y][uyx.x] != -1) road_cnt |= 1;
			}

			for (const auto& dyx : { Point(1, 0), Point(-1, 0) }) {
				auto uyx = p + dyx;
				if (board[uyx.y][uyx.x] != -1) road_cnt |= 2;
			}
			if (road_cnt == 3) {
				point_to_crossing_idx[y][x] = crossings.size();
				crossings.push(p);
			}
		}
	}

	// 交差点から見える地点を列挙
	rep(idx_crossings, crossings.size()) {
		const auto& crossing = crossings[idx_crossings];
		auto& vis = visibles[idx_crossings];
		ASSERT(vis.size() == 0, "初期化がおかしい");
		vis.push(crossing);
		for (const auto& dyx : { Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1) }) {
			auto p = crossing;
			while (true) {
				p += dyx;
				if (board[p.y][p.x] == -1) break;
				vis.push(p);
			}
		}
	}

	// 交差点同士の距離を計算
	rep(start_idx_crossing, crossings.size()) {
		Dijkstra(start_idx_crossing);
		rep(goal_idx_crossing, crossings.size()) {
			distances[start_idx_crossing][goal_idx_crossing] = dijkstra::crossing_dists[goal_idx_crossing];
		}
	}

	// 最適化
	// 可能な限り繰り返す
	{
		auto crossing_order = Stack<short, MAX_N_CROSSINGS>(crossings.size());
		iota(crossing_order.begin(), crossing_order.end(), 0);
		auto covering_set = Stack<short, MAX_N_CROSSINGS + 1>();
		auto looked = array<array<bool, 71>, 71>();

		auto best_score = (int)1e9;
		auto best_covering_set = Stack<short, MAX_N_CROSSINGS + 1>();

		auto iteration_ordering = 0;

		while (time() - t0 < 2.9) {
			iteration_ordering++;
			// 被覆集合を作るための順列を生成する (0 番目は 0 で固定)
			rep(i, crossings.size() - 2) {
				const auto r = rng.randint(1, crossings.size() - i);
				swap(crossing_order[r], crossing_order[crossings.size() - i - 1]);
			}
			covering_set.resize(0);

			// 被覆集合を作る
			auto look_cnt = 0;
			for (auto&& l : looked)for (auto&& lk : l) lk = false;
			rep(idx_order, crossing_order.size()) {
				const auto& idx_crossings = crossing_order[idx_order];
				const auto& vis = visibles[idx_crossings];
				auto f = false;  // 有用な交差点かどうか
				for (const auto& p : vis) {
					if (!looked[p.y][p.x]) {
						looked[p.y][p.x] = true;
						look_cnt++;
						f = true;
					}
				}
				if (f) covering_set.push(idx_crossings);
				if (look_cnt == n_roads) break;
			}


			// 巡回セールスマン

			auto n_nodes = covering_set.size();
			covering_set.push(0);  // 終点
			ASSERT(covering_set[0] == 0, "最初の場所は 0 のはずだよ");
			auto calc_path_dist = [&]() {
				auto res = 0;
				rep(i, n_nodes) {
					res += distances[covering_set[i]][covering_set[i + 1]];
				}
				return res;
			};

			auto n_iteration = N * N * 4;
			auto score = calc_path_dist();
			auto tmp_best_score = best_score;
			rep(iteration, n_iteration) {
				auto l = rng.randint(1, n_nodes);
				auto r = rng.randint(1, n_nodes);  // 閉区間
				if (l == r) continue;
				if (l > r) swap(l, r);
				const auto& a = covering_set[l - 1];
				const auto& b = covering_set[l];
				const auto& c = covering_set[r];
				const auto& d = covering_set[r + 1];
				const auto old_dist = (int)distances[a][b] + (int)distances[c][d];
				const auto new_dist = (int)distances[a][c] + (int)distances[b][d];
				const auto& pb = crossings[b];
				const auto& pc = crossings[c];
				const auto delta = new_dist - old_dist + board[pb.y][pb.x] - board[pc.y][pc.x];
				if (delta < 0) {
					reverse(&covering_set[l], &covering_set[r + 1]);
					ASSERT(calc_path_dist() - score == delta, "差分計算失敗！！！！");  // NDEBUG!!!!!!!!!!!!!!!!!!
					score += delta;
				}
				if (iteration_ordering > 5000 && n_iteration - 1 == iteration && tmp_best_score > score) {
					tmp_best_score = score;
					n_iteration += N * N * 4;
				}
			}

			if (best_score > score) {
				best_score = score;
				best_covering_set = covering_set;
			}
		}

		// 復元
		auto path = Stack<Point, 71 * 71 * 2>();
		path.emplace(sy, sx);
		for (int idx_covering_set = best_covering_set.size() - 2; idx_covering_set >= 0; idx_covering_set--) {
			const auto& l = best_covering_set[idx_covering_set];
			const auto& r = best_covering_set[idx_covering_set + 1];
			const auto& pr = crossings[r];
			Dijkstra(l);
			auto p = dijkstra::from[pr.y][pr.x];
			while (p != Point(-1, -1)) {
				path.push(p);
				p = dijkstra::from[p.y][p.x];
			}
		}
		reverse(path.begin(), path.end());

		// 向きに変換して出力
		rep(i, path.size() - 1) {
			const auto& pl = path[i];
			const auto& pr = path[i + 1];
			const auto d = pr - pl;
			if (d.y > 0) cout << "D";
			else if (d.y < 0) cout << "U";
			else if (d.x > 0) cout << "R";
			else cout << "L";
		}
		cout << endl;
		//cout << iteration_ordering << endl;
	}


}

int main() {
	Solve();
}
