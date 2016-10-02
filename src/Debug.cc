#include "Debug.hh"
#include "Types.hh"
#include <iostream>
#include <cstdlib> // std::abort
#include <cmath>

//http://patrick.wagstrom.net/old/misc/dotfiles/dotbash_profile.html
const char * BLACK = "\033[0;30m";
const char * RED = "\033[0;31m";
const char * GREEN = "\033[0;32m";
const char * BROWN = "\033[0;33m";
const char * BLUE = "\033[0;34m";
const char * MAGENTA = "\033[0;35m";
const char * CYAN = "\033[0;36m";
const char * WHITE = "\033[0;37m";
const char * LIGHTBLACK = "\033[1;30m";
const char * LIGHTRED = "\033[1;31m";
const char * LIGHTGREEN = "\033[1;32m";
const char * YELLOW = "\033[1;33m";
const char * LIGHTBLUE = "\033[1;34m";
const char * LIGHTMAGENTA = "\033[1;35m";
const char * LIGHTCYAN = "\033[1;36m";
const char * LIGHTWHITE = "\033[1;37m";

namespace internal {

void checkFct( bool b, const char * const expression, const std::string & message, const char * const filename, int line )
{
   if( !b )
   {
      std::cerr << "Check failed!\n\tFile:       " << filename << ":" << line << "\n" << "\tExpression: " << expression << std::endl;
      if( message.size() > 0 ) std::cerr << "\tMessage: \n\t" << message << std::endl;
   }
   std::abort();
}

void assertFct( bool b, const char * const expression, const std::string & message, const char * const filename, int line )
{
   if( !b )
   {
      std::cerr << "Assertion failed!\n\tFile:       " << filename << ":" << line << "\n" << "\tExpression: " << expression << std::endl;

      if( message.size() > 0 ) std::cerr << "\tMessage: \n\t" << message << std::endl;
   }
   std::abort();
}

void warnFct( const std::string & message, const char * const filename, int line )
{
   std::cerr << "[WARNING]" << message << std::endl << "\tFile:       " << filename << ":" << line << "\n";
}

void errorFct( const std::string & message, const char * const filename, int line )
{
   std::cerr << "[ERROR]" << message << std::endl << "\tFile:       " << filename << ":" << line << "\n";
}

void progressFct( const std::string & message )
{
   std::cout << "[PROGRESS] " << message << std::endl;
}

void checkRealEqual( real x, real y, const double epsilon, const char * const filename, int line )
{
   real absDiff = std::fabs( static_cast< long double >( x - y ) );
   if( ! ( absDiff < epsilon ) )
   {
      std::cerr << "Check for real equality failed!\n\tFile:       " << filename << ":" << line << "\n" << "\tExpression: |(" << x << ") - (" << y << ")| < "
                << epsilon << " | Result: " << absDiff << std::endl;
      std::abort();
   }
}

} // namespace internal
