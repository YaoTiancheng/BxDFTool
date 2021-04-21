#include "PCH.h"
#include "MPDispatcher.h"
#include "MPWorkerThread.h"

static const auto s_WorkerThreadCount = std::thread::hardware_concurrency();

void MP::Dispatch( uint32_t threadSize, uint32_t threadCount, LaneFunctionType laneFunction )
{
    SDispatchContext context( threadSize, threadCount, laneFunction ) ;
    context.m_NextThreadIndex = 0;

    std::vector<MP::CWorkerThread*> workerThreads;
    workerThreads.resize( s_WorkerThreadCount );
    for ( auto& wt : workerThreads )
    {
        wt = new MP::CWorkerThread( &context );
    }

    for ( auto& wt : workerThreads )
    {
        wt->Join();
        delete wt;
    }
}
