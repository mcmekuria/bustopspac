//: C05:StringConv.h
// Function templates to convert to and from strings.
#ifndef STRINGCONV_H
#define STRINGCONV_H
#include <string>
#include <sstream>
#include <cctype> 

 
template<typename T> T fromString(const std::string& s) {
  std::istringstream is(s);
  T t;
  is >> t;
  return t;
}
 
template<typename T> std::string toString(const T& t) {
  std::ostringstream s;
  s << t;
  return s.str();
}

template<typename T> T stringUpper(const std::string& s) {
  std::istringstream is(s);
  T t;
  is >> uppercase >> t;
  return t;
}

template<typename T> T stringNoUpper(const std::string& s) {
  std::istringstream is(s);
  T t;
  is >> nouppercase >> t;
  return t;
}

struct to_upper {
  int operator() ( int ch )
  {
    return std::toupper ( ch );
  }
};

struct to_lower {
  int operator() ( int ch )
  {
    return std::tolower ( ch );
  }
};
/* C++11 not working
bool is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(), 
		s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
};
*/

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
};

#endif // STRINGCONV_H ///:~
