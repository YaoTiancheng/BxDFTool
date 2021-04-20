#pragma once

namespace MP
{
    using LaneFunctionType = std::function<void( uint32_t, uint32_t, uint32_t )>;


    struct SDispatchContext
    {
        uint32_t         m_ThreadSize;
        LaneFunctionType m_LaneFunction;
    };


    class CWorkerThread
    {
    public:
        explicit CWorkerThread( const SDispatchContext* context );

        ~CWorkerThread();

        void     Deactive() { m_IsDeactived.store( true ); }

        bool     IsDeactived() { return m_IsDeactived.load(); }

        void     Join();

        void     KickOff() { m_IsDone.store( false ); }

        bool     IsDone() { return m_IsDone.load(); }

        void     SetThreadIndex( uint32_t value ) { m_ThreadIndex = value; }


    private:
        void Execute();

        const SDispatchContext* m_Context;
        std::thread*            m_Thread = nullptr;

        uint32_t                m_ThreadIndex = 0;
        std::atomic<bool>       m_IsDone;
        std::atomic<bool>       m_IsDeactived;
    };
}