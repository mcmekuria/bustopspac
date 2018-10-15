//: C06:PrintSeq.h
// Prints the contents of any sequence.
#ifndef PRINTSEQUENCE_H
#define PRINTSEQUENCE_H
#include <algorithm>
#include <iostream>
#include <iterator>
 
template<typename Iter>
void print1(Iter first, Iter last, const char* nm = "",
           const char* sep = "\n",
           std::ostream& os = std::cout) {
  if(nm != 0 && *nm != '\0')
    os << nm << " - " << sep;
  typedef typename
	  std::iterator_traits<Iter>::value_type T;
  std::copy(first, last,
            std::ostream_iterator<T>(os, sep));
  os << std::endl;
}

template<typename Iter>
void print2(Iter first, Iter last, const char* nm = "",
           const char* sep = "\n",
           std::ostream& os = std::cout) {
  if(nm != 0 && *nm != '\0')
    os << nm << " - " << sep;
  typedef typename
	  std::iterator_traits<Iter>::value_type T;
  std::copy(first, last,std::ostream_iterator<T>(os));
  os << std::endl;
}

#endif // PRINTSEQUENCE_H ///:~

// Prints a sequence container
#ifndef PRINT_CONTAINER_H
#define PRINT_CONTAINER_H
 
template<class Cont>
void print(Cont& c, const char* nm = "",
           const char* sep = "\n",
           std::ostream& os = std::cout) {
  print(c.begin(), c.end(), nm, sep, os);
}
#endif ///:~
