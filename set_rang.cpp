#include <set>
#include <iostream>
#include "linkedlist.h"

using namespace std;

typedef set<double, less<double> > SET_Book;

void main(void)
 {
   linked_list<Book> list;
   Book cbib("Jamsa's C/C++/C# Programmer's Bible", "Jamsa", "Delmar", 49.95);
   Book vbtips("1001 Visual Basic Programmer's Tips", "Jamsa and Klander", "Jamsa Press", 54.95);
   Book hacker("Hacker Proof", "Klander", "Jamsa Press", 44.95);
   Book c;

   list_object<Book> *p;

   list.store(cbib);
   list.store(vbtips);
   list.store(hacker);
  
 
   SET_Book s1;
   SET_Book::iterator i;


   s1.insert(cbib.get_price());
   s1.insert(vbtips.get_price());
   s1.insert(hacker.get_price());
   cout << "s1 -- starting at s1.lower_bound(45)" << endl<< endl;
   for (i = s1.lower_bound(45); i != s1.end(); i++)
      cout << "s1 has " << *i << " in its set." << endl;
     cout <<endl;
   cout << "s1 -- starting at s1.lower_bound(50)" << endl;
// prints: 15,20,25
   for (i = s1.lower_bound(50);i != s1.end(); i++)
      cout << "s1 has " << *i << " in its set." << endl;
     cout <<endl;

   cout << "s1 -- starting at s1.upper_bound(55)" << endl;
// prints: 15,20,25
   for (i = s1.upper_bound(55); i != s1.end(); i++)
      cout << "s1 has " << *i << " in its set." << endl;
     cout <<endl;

	 cout << "s1 -- starting at s1.upper_bound(55)" << endl;
// prints: 20,25
   for (i = s1.upper_bound(50); i != s1.end(); i++)
     cout << "s1 has " << *i << " in its set." << endl;
     cout <<endl;

   cout << "s1 -- s1.equal_range(45)" << endl;
// does not print anything
   for (i = s1.equal_range(45).first;i != s1.equal_range(45).second; i++)
      cout << "s1 has " << *i << " in its set." << endl;
     cout <<endl;

   cout << "s1 -- s1.equal_range(50)" << endl;
// prints: 15
   for (i = s1.equal_range(50).first;i != s1.equal_range(50).second; i++)
      cout << "s1 has " << *i << " in its set." << endl;
     cout <<endl;


   cout << "'Manually' walk through the list." << endl;
   p = list.getstart();
   while(p) {
      p->getinfo(c);
      c.show();
      p = p->getnext();
    }
   cout << endl << endl;
 
}
