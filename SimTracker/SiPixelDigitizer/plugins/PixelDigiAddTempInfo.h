#ifndef PixelDigiAddTempInfo_h
#define PixelDigiAddTempInfo_h

#include "DataFormats/GeometryVector/interface/LocalPoint.h"

class PixelDigiAddTempInfo {
public:
  PixelDigiAddTempInfo(
      unsigned int ch, size_t Hindex, Local3DPoint entryP , Local3DPoint exitP) {
    chan = ch;
    index = Hindex;
    TheEntryPoint = entryP;
    TheExitPoint = exitP;
  };
  PixelDigiAddTempInfo() {
    chan = 0;
    index = 0;
//    TheEntryPoint(0,0,0);
//    TheExitPoint(0,0,0);
  };
  ~PixelDigiAddTempInfo(){};
  unsigned int channel() const { return chan; };
  size_t hitIndex() const { return index; };
  Local3DPoint entryPoint() const { return TheEntryPoint; };
  Local3DPoint exitPoint() const { return TheExitPoint; }

  inline bool operator<(const PixelDigiAddTempInfo& other) const { return channel() < other.channel(); }

private:
  unsigned int chan;
  size_t index;
  Local3DPoint TheEntryPoint;
  Local3DPoint TheExitPoint;
};
#endif

