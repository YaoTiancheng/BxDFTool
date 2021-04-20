#include "PCH.h"
#include "MPDispatcher.h"
#include "MPWorkerThread.h"

static const auto s_WorkerThreadCount = std::thread::hardware_concurrency();

void MP::Dispatch( uint32_t threadSize, uint32_t threadCount, LaneFunctionType laneFunction )
{
    std::list<uint32_t> remainingThreadIndices;
    for ( uint32_t i = 0; i < threadCount; ++i )
    {
        remainingThreadIndices.push_back( i );
    }

    SDispatchContext context;
    context.m_ThreadSize   = threadSize;
    context.m_LaneFunction = laneFunction;

    // Start all the worker threads
    std::vector<MP::CWorkerThread*> workerThreads;
    workerThreads.resize( s_WorkerThreadCount );
    for ( auto& wt : workerThreads )
    {
        wt = new MP::CWorkerThread( &context );
    }

    while ( !remainingThreadIndices.empty() )
    {
        for ( auto& wt : workerThreads )
        {
            if ( !wt->IsDeactived() && wt->IsDone() )
            {
                // Kickoff a new thread task if available
                if ( !remainingThreadIndices.empty() )
                {
                    uint32_t threadIndex = remainingThreadIndices.front();
                    wt->SetThreadIndex( threadIndex );
                    wt->KickOff();

                    remainingThreadIndices.erase( remainingThreadIndices.begin() );
                }
                else
                {
                    wt->Deactive();
                }
            }
        }
    }

    for ( auto& wt : workerThreads )
    {
        wt->Join();
        delete wt;
    }
}
