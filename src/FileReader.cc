#include <string>
#include <algorithm>
#include "FileReader.hh"

FileReader::FileReader()
{
    reset();
}

void FileReader::reset()
{
    realParams_.clear();
    intParams_.clear();
    strParams_.clear();
    const std::map< std::string, real > realParameters = {
        /* Domain size */
        { "xlength", real(1) }, { "ylength", real(1) },
        /* Reynolds number */
        { "re", real(1) },
        /* Timestep size and safety factor (tao) */
        { "dt", real(0.01) }, { "safetyfactor", real(1) },
        /* External forces */
        { "gx", real(0) }, { "gy", real(0) },
        /* Velocities at the boundaries */
        { "boundary_velocity_n", real(0) }, { "boundary_velocity_s", real(0) },
        { "boundary_velocity_w", real(0) }, { "boundary_velocity_e", real(0) },
        /* Initial values for U, V, P */
        { "u_init", real(0) }, { "v_init", real(0) }, { "p_init", real(0) },
        /* Epsilon, Omega (SOR relaxation parameter), Gamma */
        { "eps", real(0) }, { "omg", real(0) }, { "gamma", real(0.1) },
        /* Rectangular obstacle */
        { "rectanglex1", real(0) }, { "rectangley1", real(0) }, { "rectanglex2", real(0) }, { "rectangley2", real(0) },
        /* Circle obstacle */
        { "circlex", real(0) }, { "circley", real(0) }, { "circler", real(0) }
    };

    const std::map< std::string, int > intParameters = {
        /* Domain mesh size in x, y directions */
        { "imax", int(1)  }, { "jmax", int(1) },
        /* Nr of timesteps */
        { "timesteps", int(100) },
        /* Boundary conditions for all directions */
        { "boundary_condition_n", NOSLIP }, { "boundary_condition_s", NOSLIP },
        { "boundary_condition_w", NOSLIP }, { "boundary_condition_e", NOSLIP },
        /* Number of iterations */
        { "itermax", int(1) },
        /* Residual check frequency */
        { "checkfrequency", int(1) },
        /* Pressure normalization frequency */
        { "normalizationfrequency", int(1) },
        /* VTK output frequency */
        { "outputinterval", int(1) },
        /* Flag to describe the type of obstacle */
        { "obstacletype", int(NO_OBSTACLE) }
    };

    const std::map< std::string, std::string > strParameters = {
        { "name", "" }
    };

    for( auto it = realParameters.begin(); it != realParameters.end(); ++it )
        registerRealParameter( it->first, it->second );

    for( auto it = intParameters.begin(); it != intParameters.end(); ++it )
        registerIntParameter( it->first, it->second );

    for( auto it = strParameters.begin(); it != strParameters.end(); ++it )
        registerStringParameter( it->first, it->second );
}

void FileReader::registerIntParameter(const std::string &key, int init)
{
    intParams_.insert( std::make_pair( key , init) );
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
    realParams_.insert( std::make_pair( key , init) );
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
    strParams_.insert( std::make_pair( key , init) );
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
    auto search = strParams_.find(key);
    if ( search != strParams_.end() )
    {
        search->second = in;
    }
    else
    {
        WARN( "Parameter '" + key + "' not found, creating new.");
        strParams_.insert( std::make_pair( key, in ) );
    }
}

void FileReader::setParameter(const std::string &key, real in)
{
    auto search = realParams_.find(key);
    if ( search != realParams_.end() )
    {
        search->second = in;
    }
    else
    {
        WARN( "Parameter '" + key + "' not found, creating new.");
        realParams_.insert( std::make_pair( key, in ) );
    }
}

void FileReader::setParameter(const std::string &key, int in)
{
    auto search = intParams_.find(key);
    if ( search != intParams_.end() )
    {
        search->second = in;
    }
    else
    {
        WARN( "Parameter '" + key + "' not found, creating new.");
        intParams_.insert( std::make_pair( key, in ) );
    }
}


bool FileReader::readFile(const std::string &name)
{
    CHECK( !name.empty() );
    std::ifstream paramFile( name );
    std::string line;

    if ( !paramFile )
    {
        std::cerr << "Failed to open the configuration file!" << std::endl;
        return false;
    }

    // reset simulation parameters before reading the new file
    reset();

    while( getline(paramFile, line) )
    {
        if ( line.empty() )
        {
            continue;
        }

        std::string word, paramName;
        std::istringstream lineStream(line);
        bool isParamValue = false;

        while( lineStream >> word )
        {
            if ( word[0] == '#' )
            {
                break;
            }

            if (isParamValue)
            {
                auto intSearch = intParams_.find( paramName );
                if ( intSearch != intParams_.end() )
                {
                    int intParam;
                    if ( word == "noslip" )
                        intParam = NOSLIP;
                    else if ( word == "slip" )
                        intParam = SLIP;
                    else if ( word == "outflow" )
                        intParam = OUTFLOW;
                    else if ( word == "inflow" )
                        intParam = INFLOW;
                    else if ( word == "periodic" )
                        intParam = PERIODIC;
                    else
                        intParam = std::stoi( word );

                    setParameter( paramName, intParam );
                    isParamValue = false;
                    continue;
                }

                auto realSearch = realParams_.find( paramName );
                if ( realSearch != realParams_.end() )
                {
                    real realParam = static_cast<real>( std::stod( word ) );
                    setParameter( paramName, realParam );
                    isParamValue = false;

                    // check for obstacle
                    if ( (paramName == "rectanglex1") || (paramName == "rectangley1") || (paramName == "rectanglex2") || (paramName == "rectangley2") ) {
                        setParameter( "obstacletype", int(RECTANGLE) );
                    } else if ( (paramName == "circlex") || (paramName == "circley") || (paramName == "circler") ) {
                        setParameter( "obstacletype", int(CIRCLE) );
                    }
                    continue;
                }

                setParameter( paramName, word );
                isParamValue = false;
            }
            else
            {
                std::transform(word.begin(), word.end(), word.begin(), ::tolower); // transform string to lower case
                paramName = word;
                isParamValue = true;
            }
        }
    }

    return true;
}



void FileReader::printParameters() const
{
    std::cout << "Simulation Parameters" << std::endl;
    for( auto it = realParams_.begin(); it != realParams_.end(); ++it )
        std::cout << it->first << " = " << it->second << std::endl;

    for( auto it = intParams_.begin(); it != intParams_.end(); ++it )
        std::cout << it->first << " = " << it->second << std::endl;

    for( auto it = strParams_.begin(); it != strParams_.end(); ++it )
        std::cout << it->first << " = " << it->second << std::endl;
}

inline bool FileReader::isInteger(const std::string & s)
{
    if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+') && ( s.find(".") != std::string::npos ) )) return false;

    char *p;
    strtol( s.c_str(), &p, 10 ) ;

    return (*p == 0);
}

inline bool FileReader::isReal(const std::string & s)
{
    if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false;

    char *p;
    strtod( s.c_str(), &p );

    return (*p == 0);
}
