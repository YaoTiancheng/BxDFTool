#pragma once

namespace MP
{
    using LaneFunctionType = std::function<void( uint32_t, uint32_t, uint32_t )>;


    struct SDispatchContext
    {
        SDispatchContext( uint32_t threadSize, uint32_t threadCount, LaneFunctionType laneFunction ) 
            : m_ThreadSize( threadSize ), m_ThreadCount( threadCount ), m_LaneFunction( laneFunction ), m_NextThreadIndex( 0 )
        {}

        const uint32_t          m_ThreadSize;
        const uint32_t          m_ThreadCount;
        const LaneFunctionType  m_LaneFunction;

        std::atomic_uint32_t    m_NextThreadIndex;
    };


    class CWorkerThread
    {
    public:
        explicit CWorkerThread( SDispatchContext* context );

        ~CWorkerThread();

        void     Join();

    private:
        void     Execute();

        SDispatchContext*       m_Context;
        std::thread*            m_Thread = nullptr;
    };
}