
#ifndef _BSGTREE_HPP_
#define _BSGTREE_HPP_

#include "quadrics.hpp"

typedef Geometry::Quadric<2, float> q2f;
typedef Geometry::Quadric<3, float> q3f;

typedef Geometry::Quadric<2, double> q2d;
typedef Geometry::Quadric<3, double> q3d;

class BSGNode {
public:
  BSGNode() {
    parent = NULL;
    refs = 0;
  }
  
  BSGNode(const BSGNode &src) {
    refs = 0;
    addRef(parent, src.parent);
    addRef(sibilings[0], src.sibilings[0]);
    addRef(sibilings[1], src.sibilings[1]);
    addRef(children, src.children);
  }
  
  virtual ~BSGNode() {
    remRef(parent);
    remRef(sibilings[0]);
    remRef(sibilings[1]);
    remRef(children);
  }
  
  bool isRoot() {
    if(parent == NULL)
      return true;
    else
      return false;
  }
protected:
  q3f quadric;
  unsigned refs;
  BSGNode *parent;
  BSGNode *sibilings[2];
  BSGNode *children;
  
  static void addRef(BSGNode * &newRef, BSGNode *ref) {
    newRef = ref;
    if(newRef)
      newRef->refs++;
  }
  
  static void remRef(BSGNode * &oldRef) {
    if(oldRef) {
      oldRef->refs--;
      if(oldRef->refs <= 0)
        delete oldRef;
      oldRef = NULL;
    }
  }
};

#endif
