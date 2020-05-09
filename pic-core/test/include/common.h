#ifndef COMMON_H
#define COMMON_H

#include <catch2.hpp>
#include <cstdlib>
#include <type_traits>
#include <typeinfo>
#include <cstdlib>

#include <csignal>
#include <iostream>
#include <ctime>
#include <chrono>

template <class T>
constexpr
std::string_view
type_name()
{
    using namespace std;
#ifdef __clang__
    string_view p = __PRETTY_FUNCTION__;
    return string_view(p.data() + 34, p.size() - 34 - 1);
#elif defined(__GNUC__)
    string_view p = __PRETTY_FUNCTION__;
#  if __cplusplus < 201402
    return string_view(p.data() + 36, p.size() - 36 - 1);
#  else
    return string_view(p.data() + 49, p.find(';', 49) - 49);
#  endif
#elif defined(_MSC_VER)
    string_view p = __FUNCSIG__;
    return string_view(p.data() + 84, p.size() - 84 - 7);
#endif
}

#define TYPE_NAME(type) type_name<decltype(type)>()
using namespace Catch::literals;

#define REQUIRE3d(_a, _b){\
auto a = _a;\
auto b = _b;\
REQUIRE(Approx(a[0]) == b[0]);\
REQUIRE(Approx(a[1]) == b[1]);\
REQUIRE(Approx(a[2]) == b[2]);\
}

//Time Funcs
using ns=std::chrono::nanoseconds;
using tp=std::chrono::high_resolution_clock::time_point;
constexpr auto tnow=std::chrono::high_resolution_clock::now;

ns elapsed(tp start){
	return std::chrono::duration_cast<ns> (tnow() - start);
}

std::ostream &operator<<(std::ostream &out, const std::chrono::nanoseconds & t){
    out<<t.count()<<"[ns]";
    return out;
}


auto timeFuncInvocation = 
    [](auto&& func, auto&&... params) {
        const auto& start = tnow();
        std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
        return elapsed(start);
     };

static volatile sig_atomic_t sig_caught = 0;

void handle_sighup(int signum) 
{
    /* in case we registered this handler for multiple signals */ 
    if (signum == SIGHUP) {
        sig_caught = 1;
    }
    else if (signum == SIGABRT) {
        sig_caught = 2;
    }
}

#endif /* COMMON_H */