#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>

#include <SDL_timer.h>

#include "interface/window.h"
#include "macros/finally.h"
#include "utils/clock.h"
#include "utils/metronome.h"
#include "program/platform.h"

#if IMP_PLATFORM_IS(emscripten)
#include <emscripten.h>
#endif

namespace Program
{
    class BasicState
    {
        bool loop_running = false;

      public:
        virtual void BeginFrame() {}
        virtual void EndFrame() {}

        virtual void Tick() {}
        virtual void Render() {}

        // Should call `BeginFrame()` once, then `Tick()` and `Frame()` 0 or more times, then `EndFrame()`.
        // Should return `false` to indicate that the loop should be stopped.
        virtual bool RunSingleFrame() = 0;

        void RunMainLoop()
        {
            // Protect against double entry.
            if (loop_running)
                return;
            loop_running = true;
            FINALLY{loop_running = false;};

            // Run the loop.
            #if !IMP_PLATFORM_IS(emscripten)
            while (RunSingleFrame()) {}
            #else
            emscripten_set_main_loop_arg([](void *self){static_cast<BasicState *>(self)->RunSingleFrame();}, this, -1, true);
            #endif
        }
    };

    class DefaultBasicState : public BasicState
    {
        bool executing_frame = false;
        std::uint64_t frame_start = -1;

      protected:
        bool stop = false;

      public:
        // Returns the tick metronome, or `nullptr` if want to run one tick per frame.
        // The intended way of overriding is to add a metronome as a class member, and return a pointer to it.
        virtual Metronome *GetTickMetronome() {return nullptr;}

        // Returns true if it's advisable to have a FPS cap. That is, when vsync is disabled.
        // See comment on `GetFpsCap` for the intended use of this function.
        [[nodiscard]] static bool NeedFpsCap()
        {
            return Interface::Window::IsOpen() && Interface::Window::Get().VSyncMode() == Interface::VSync::disabled;
        }

        // Returns target FPS. `<= 0` if not limited.
        // When overriding this function, it's recommended to multiply the return value by `NeedFpsCap()`,
        // which enables the cap only when vsync is disabled.
        virtual int GetFpsCap() {return 0;}

        // Ignored if FPS cap is disabled.
        // FPS is capped by adding a delay after frames that are too short.
        // The delay will be partially created using a `sleep` function, and partially using a busy loop.
        // This function returns the preferred duration of the busy loop.
        // If it returns -1, the busy loop will not be used. It lowers CPU load, but might make the timing less precise.
        // If it returns 0, the busy loop will be used as little as possible.
        // If it returns a large value, `sleep` will not be used at all, which increases CPU loads but improves precision.
        virtual int GetFpsCapPreferredBusyLoopDurationMs() {return 1;}

        bool RunSingleFrame() override
        {
            if (executing_frame)
                return !stop;
            executing_frame = true;
            FINALLY{executing_frame = false;};

            // Load some basic config from state.
            auto *metronome = GetTickMetronome();
            auto fps_cap = GetFpsCap();
            bool have_fps_cap = fps_cap > 0 && !IMP_PLATFORM_IS(emscripten);

            // Compute timings if needed.
            std::uint64_t delta = 0;
            if (metronome || have_fps_cap)
            {
                std::uint64_t new_frame_start = Clock::Time();

                if (frame_start != std::uint64_t(-1))
                    delta = new_frame_start - frame_start;

                frame_start = new_frame_start;
            }

            // Begin frame.
            BeginFrame();

            // Tick.
            if (metronome)
            {
                while (metronome->Tick(delta))
                    Tick();
            }
            else
            {
                Tick();
            }

            // Render.
            Render();

            // End frame.
            EndFrame();

            // Cap FPS.
            if (have_fps_cap)
            {
                std::uint64_t clock_ticks_per_sec = Clock::TicksPerSecond();
                std::uint64_t desired_frame_len = 1.0 / fps_cap * clock_ticks_per_sec;
                std::uint64_t frame_end = Clock::Time(), frame_len = frame_end - frame_start;

                // If necessary, add a delay.
                if (frame_len < desired_frame_len)
                {
                    int sleep_ms = (desired_frame_len - frame_len) * 1000 / clock_ticks_per_sec;

                    auto busy_loop_len_ms = GetFpsCapPreferredBusyLoopDurationMs();
                    bool busy_loop_allowed = busy_loop_len_ms >= 0;

                    if (busy_loop_allowed)
                        sleep_ms -= busy_loop_len_ms;

                    // Sleep.
                    if (sleep_ms > 0)
                        SDL_Delay(sleep_ms);

                    // Busy loop.
                    if (busy_loop_allowed)
                    {
                        std::uint64_t busy_loop_end = frame_start + desired_frame_len;
                        while (Clock::Time() < busy_loop_end) {}
                    }
                }
            }

            return !stop;
        }
    };
}
