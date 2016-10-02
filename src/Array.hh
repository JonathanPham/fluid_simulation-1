#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include <memory>
#include "Types.hh"
#include "Debug.hh"

//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
 *
 *    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
 */
//*******************************************************************************************************************
template< typename T >
class Array
{
    public:
        // Constructors for 1D, 2D and 3D
        Array() : xLength_( 0 ), yLength_( 0 ), zLength_( 0 ), size_( 0 ), arrayData_( nullptr ) {}
        Array( int xSize ) : xLength_( xSize ), yLength_( 1 ), zLength_( 1 ), size_( xSize ) { createArrayData(); }
        Array( int xSize, int ySize ) : xLength_( xSize ), yLength_( ySize ), zLength_( 1 ), size_( xSize * ySize ) { createArrayData(); }
        Array( int xSize, int ySize, int zSize ) : xLength_( xSize ), yLength_( ySize ), zLength_( zSize ), size_( xSize * ySize * zSize ) { createArrayData(); }
        // copy-constructor
        Array( const Array & s ) : xLength_( s.getSize( 0 ) ), yLength_( s.getSize( 1 ) ), zLength_( s.getSize( 2 ) ), size_( s.getSize() )
        {
            createArrayData();
            for( int i = 0; i < s.getSize(); ++i ) operator()(i) = s(i);
        }

        Array operator=(const Array & s);

        // Access Operators for 1D, 2D and 3D
        inline T & operator ()( int i );
        inline T & operator ()( int i, int j );
        inline T & operator ()( int i, int j, int k );
        // for const Arrays the following access operators are required
        inline const T & operator ()( int i ) const;
        inline const T & operator ()( int i, int j ) const;
        inline const T & operator ()( int i, int j, int k ) const;

        // initialize the whole array with a constant value
        void fill( T value );

        // return total size of the array
        int getSize() const { return size_; }

        // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
        // other dimension values are not allowed
        int getSize( int dimension ) const;

        // Print the whole array ( for debugging purposes )
        void print();

    private:
        int xLength_;
        int yLength_;
        int zLength_;
        int size_;
        std::unique_ptr< T[] > arrayData_; // smart pointer to array data

        inline void createArrayData();
};


// copy-assignment operator
template< typename T >
Array< T > Array< T >::operator=( const Array< T > & s )
{
    xLength_ = s.getSize( 0 );
    yLength_ = s.getSize( 1 );
    zLength_ = s.getSize( 2 );
    size_ = s.getSize();
    createArrayData();
    for( int i = 0; i < s.getSize(); ++i )
        operator()(i) = s( i );
    return *this;
}


// operator() 1D
template< typename T >
inline T & Array< T >::operator ()( int i )
{
    ASSERT_MSG( i < getSize(), "Access out-of-bounds! Requested index = " << i << ". Array length = " << getSize() );
    return arrayData_[i];
}


// operator() 2D
template< typename T >
inline T & Array< T >::operator ()( int i, int j )
{
    const int index = j * xLength_ + i;
    ASSERT_MSG( index < getSize(), "Access out-of-bounds! Requested index = " << index << " (i = " << i << ", j = " << j << "). Array length = " << getSize() );
    return arrayData_[index];
}


// operator() 3D
template< typename T >
inline T & Array< T >::operator ()( int i, int j, int k )
{
    const int index = k * xLength_ * yLength_ + j * xLength_ + i;
    ASSERT_MSG( index < getSize(), "Access out-of-bounds! Requested index = " << index << ". Array length = " << getSize() );
    return arrayData_[index];
}


// const operator() 1D
template< typename T >
inline const T & Array< T >::operator ()( int i ) const
{
    ASSERT_MSG( i < getSize(), "Access out-of-bounds! Requested index = " << i << ". Array length = " << getSize() );
    return arrayData_[i];
}


// const operator() 2D
template< typename T >
inline const T & Array< T >::operator ()( int i, int j ) const
{
    const int index = j * xLength_ + i;
    ASSERT_MSG( index < getSize(), "Access out-of-bounds! Requested index = " << index << ". Array length = " << getSize() );
    return arrayData_[index];
}


// const operator() 3D
template< typename T >
inline const T & Array< T >::operator ()( int i, int j, int k ) const
{
    const int index = k * xLength_ * yLength_ + j * xLength_ + i;
    ASSERT_MSG( index < getSize(), "Access out-of-bounds! Requested index = " << index << ". Array length = " << getSize() );
    return arrayData_[index];
}


// allocate array data
template< typename T >
inline void Array< T >::createArrayData()
{
    arrayData_ = std::unique_ptr< T[] >( new T[reinterpret_cast< int >( getSize() )] );
}


// Print the whole array (for debugging purposes)
template< typename T >
void Array< T >::print()
{
    // For 2D Arrays the positive x-coordinate goes to the right
    //                   positive y-coordinate goes upwards
    //      -> the line with highest y-value should be printed first
    for( int z = 0; z < zLength_; ++z )
    {
        std::cout << "z = " << z << std::endl;
        for( int y = yLength_ - 1; y >= 0; --y )
        {
            for( int x = 0; x < xLength_; ++x )
                std::cout << operator()( x, y, z ) << " ";
            std::cout << std::endl;
        }
    }
}


template< typename T >
void Array< T >::fill( T value )
{
    int dim = getSize();
    for( int i = 0; i < dim; ++i )
        operator()(i) = value;
}


template< typename T >
int Array< T >::getSize( int dimension ) const
{
    ASSERT( ( dimension >= 0 ) && ( dimension <= 2 ) );
    switch( dimension )
    {
        case 0:
            return xLength_;
            break;
        case 1:
            return yLength_;
            break;
        case 2:
            return zLength_;
            break;
        default:  // should not be reached
            return 0;
            break;
    }
}

#endif //ARRAY_HH
