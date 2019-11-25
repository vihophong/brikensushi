// erasing from map
#include <iostream>
#include <map>
#include <vector>
#include <iostream>

int testmap ()
{
  std::map<char,int> mymap;
  std::map<char,int>::iterator mymap_it;
  std::map<char,int>::iterator mymap_itp;

  // insert some values:
  mymap['a']=10;
  mymap['b']=20;
  mymap['c']=30;
  mymap['d']=40;
  mymap['e']=50;
  mymap['f']=60;

  mymap_it=mymap.find('b');
  mymap_itp=mymap_it;
  mymap_itp++;
  mymap.erase (mymap.begin(),mymap_itp);                   // erasing by iterator
  cout<<mymap.size()<<endl;
  // show content:
  //std::cout <<  mymap_it->first << " => " <<  mymap_it->second << '\n';

  return 0;
}
