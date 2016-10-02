#ifndef DEBUG_HH
#define DEBUG_HH

#include <sstream>
#include <limits>
#include "Types.hh"

//===================================================================================================================
//
//  CHECK Macro used for tests, is activated in Debug and Release Mode
//
//===================================================================================================================

#define CHECK_MSG(X, MSG) \
   if( !(X) )  \
   { std::stringstream ss; \
     ss << MSG; \
     internal::checkFct ( (X), #X, ss.str(), __FILE__, __LINE__ );\
   }

#define CHECK(X) \
   if( !(X) ) {  internal::checkFct ( (X), #X, "", __FILE__, __LINE__ ); }

#define EPSILON (real(1e-12)) // numeric precision for comparing real numbers

#define CHECK_REAL_EQUAL(X, Y) internal::checkRealEqual( X, Y, EPSILON, __FILE__, __LINE__ );


#define PROGRESS(MSG) \
  { std::stringstream ss; \
    ss << MSG; \
    internal::progressFct( ss.str() );\
  }

#define WARN(MSG) \
   { std::stringstream ss; \
     ss << MSG; \
     internal::warnFct ( ss.str(), __FILE__, __LINE__ );\
   }

#define ERROR(MSG) \
  { std::stringstream ss; \
    ss << MSG; \
    internal::errorFct( ss.str(), __FILE__, __LINE__ );\
  }

//===================================================================================================================
//
//  ASSERT Macro checks the given expression in Debug mode, disabled in Release mode
//
//===================================================================================================================

#ifndef NDEBUG

#define ASSERT_MSG(X, MSG) \
   if( !(X) )  \
   { std::stringstream ss; \
     ss << MSG; \
     internal::assertFct ( (X), #X, ss.str(), __FILE__, __LINE__ );\
   }

#define ASSERT(X) \
   if( !(X) ) {  internal::assertFct ( (X), #X, "", __FILE__, __LINE__ ); }

#else

#define ASSERT_MSG(X, MSG)
#define ASSERT(X)

#endif //NDEBUG

namespace internal {

void checkFct( bool b, const char * const expression, const std::string & message, const char * const filename, int line );
void assertFct( bool b, const char * const expression, const std::string & message, const char * const filename, int line );

void progressFct( const std::string & message );
void warnFct( const std::string & message, const char * const filename, int line );
void errorFct( const std::string & message, const char * const filename, int line );

void checkRealEqual( real x, real y, const double epsilon, const char * const filename, int line );

}

#endif // DEBUG_HH
