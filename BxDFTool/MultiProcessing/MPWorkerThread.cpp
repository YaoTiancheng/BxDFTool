#include "PCH.h"
#include "MPWorkerThread.h"

MP::CWorkerThread::CWorkerThread( SDispatchContext* context )
    : m_Context( context )
{
    m_Thread = new std::thread( &CWorkerThread::Execute, this );
}

MP::CWorkerThread::~CWorkerThread()
{
    delete m_Thread;
}

void MP::CWorkerThread::Join()
{
    m_Thread->join();
}

void MP::CWorkerThread::Execute()
{
    while ( true )
    {
        uint32_t threadIndex = m_Context->m_NextThreadIndex.fetch_add( 1 );
        bool isQueueEmpty = ( threadIndex >= m_Context->m_ThreadCount );
        if ( isQueueEmpty )
        {
            break;
        }

        uint32_t laneCount = m_Context->m_ThreadSize;
        uint32_t laneBegin = threadIndex * laneCount;
        uint32_t laneEnd = laneBegin + laneCount;
        for ( uint32_t i = 0; i < laneCount; ++i )
        {
            m_Context->m_LaneFunction( threadIndex, i, laneBegin + i );
        }
    }
}
