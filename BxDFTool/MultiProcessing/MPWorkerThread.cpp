#include "PCH.h"
#include "MPWorkerThread.h"

MP::CWorkerThread::CWorkerThread( const SDispatchContext* context )
    : m_Context( context )
    , m_IsDone( true )
    , m_IsDeactived( false )
{
    m_Thread = new std::thread( &CWorkerThread::Execute, this );
}

MP::CWorkerThread::~CWorkerThread()
{
    delete m_Thread;
}

void MP::CWorkerThread::Join()
{
    Deactive();
    m_Thread->join();
}

void MP::CWorkerThread::Execute()
{
    while ( !m_IsDeactived.load() )
    {
        while ( m_IsDone.load() )
        {
            if ( m_IsDeactived.load() )
                return;
        }

        uint32_t laneCount = m_Context->m_ThreadSize;
        uint32_t laneBegin = m_ThreadIndex * laneCount;
        uint32_t laneEnd = laneBegin + laneCount;
        for ( uint32_t i = 0; i < laneCount; ++i )
        {
            m_Context->m_LaneFunction( m_ThreadIndex, i, laneBegin + i );
        }

        m_IsDone.store( true );
    }
}
