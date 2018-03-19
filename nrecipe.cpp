#include<cmath>
#include<limits>
#include<string>
#include<iostream>

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T SQR(const T a) {return a*a;}

inline void nrerror(const std::string s)
{ std::cerr << s << std::endl; exit(1); }
