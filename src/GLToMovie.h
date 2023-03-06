#ifndef __GLTOMOVIE_H_9AAF23E7_E751_4d77_96E2_880C6FB48B3D__
#define __GLTOMOVIE_H_9AAF23E7_E751_4d77_96E2_880C6FB48B3D__

#include "gl\gl.h"
#include "AviFile.h"

/// <Summary>
/// Creates a Movie from the OpenGL Rendering Surface Bits.
/// Movie Frame's Width and Height can not be changed once the Movie object is created.
/// Any number of objects of this class can be present simultaneously with in a single app.
/// </Summary>
class CGLToMovie
{
	CAviFile	m_MovieFile;		/*Movie File Object*/
	int			m_nWidth;			/*Movie Frame Width*/
	int			m_nHeight;			/*Movie Frame Height*/
	int			m_nBitsPerPixel;	/*Number of Bits per Pixel*/
	LPVOID		m_pBits;			/*Pointer to Memory required to hold the Frame Content*/
public:

	CGLToMovie(LPCTSTR lpszOutputMovieFileName = _T("Output.avi"), 
		int nFrameWidth = GetSystemMetrics(SM_CXSCREEN),	/*Movie Frame Width*/
		int nFrameHeight = GetSystemMetrics(SM_CYSCREEN),	/*Movie Frame Height*/
		int nBitsPerPixel = 24,		/*Bits Per Pixel*/
		DWORD dwCodec = mmioFOURCC('M','P','G','4'),	/*Video Codec for Frame Compression*/
		DWORD dwFrameRate = 1)		/*Frames Per Second (FPS)*/
		: m_MovieFile(lpszOutputMovieFileName, dwCodec, dwFrameRate)
	{
		m_nWidth = nFrameWidth;
		m_nHeight = nFrameHeight;
		m_nBitsPerPixel = nBitsPerPixel;

		m_pBits = malloc(m_nWidth * m_nHeight * m_nBitsPerPixel/8);	/*create the required memory to hold the frame content*/
	}

	~CGLToMovie(void)
	{
		if(m_pBits)
		{
			free(m_pBits);
			m_pBits = NULL;
		}
	}

	/// <Summary>
	/// Reads the content of the OpenGL Frame buffer and inserts the read content
	/// as a new frame into the movie.
	/// Should be called before SwapBuffers() for each rendered frame.
	/// </Summary>
	inline HRESULT RecordFrame()
	{
		glFlush(); glFinish();
		glReadPixels(0, 0, m_nWidth, m_nHeight, GL_BGR_EXT, GL_UNSIGNED_BYTE, m_pBits);
		return m_MovieFile.AppendNewFrame(m_nWidth, m_nHeight, m_pBits, m_nBitsPerPixel);
	}

	inline const CAviFile* MovieFile() const	{	return &m_MovieFile;	}
};

#endif