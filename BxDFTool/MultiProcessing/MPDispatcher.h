#pragma once

namespace MP
{
    using LaneFunctionType = std::function<void( uint32_t, uint32_t, uint32_t )>;

    void Dispatch( uint32_t threadSize, uint32_t threadCount, LaneFunctionType laneFunction );
}