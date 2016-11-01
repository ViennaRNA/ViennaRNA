/*
  Last changed Time-stamp: <2007-07-10 19:26:47 xtof>
  $Id: template_utils.c,v 1.2 2007/10/01 08:29:47 Kinwalker Exp $
*/

/*
 * template hocus-pocus to dump pair<T,T> to stdout using
 * algorithm copy
 */

// construct helper object
template <class T>
struct place_holder
{
  T item_held;
  place_holder(const T& item) : item_held(item) {}
};

// extract object from helper object and dump it to stream
template<class T>
std::ostream& operator<<(std::ostream& os, const place_holder<T>& t)
{
  os << t.item_held;

  return (os);
}

// dump formated pair to a stream
template<class T1, class T2>
std::ostream& operator<< (std::ostream &os, std::pair<T1,T2> p) {
  os << "[" << p.first << "," << p.second << "]";
  return os;
}

// for the following function we needed the above template hocus-pocus
template<class T1, class T2>
void
print_vector_of_pairs (std::vector<std::pair<T1,T2> >& V)
{
  copy(V.begin(), V.end(),
       std::ostream_iterator<place_holder<std::pair<T1,T2> > >(std::cout, " "));
  std::cout << std::endl;
}

/* End of file */
