#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <atomic>
#include <exception>
#include <functional>
#include <limits>
#include <ostream>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <iostream>

struct Range
{
    Range(size_t _start, size_t _end):
        start(_start),
        end(_end)
    {

    }
    size_t start;
    size_t end;

    size_t size() const
    {
        return end-start;
    }
};

class ThreadPool
{
public:
    ThreadPool(ThreadPool& other) = delete;
    ThreadPool& operator=(const ThreadPool& other) = delete;

    ~ThreadPool();

    static ThreadPool* i();
#ifndef ENQUEUE_IS_RUN
    template<typename F, typename... Args>
    void enqueue(F&& f, Args&&... args)
    {
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_tasks.emplace(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
        }

        m_cv.notify_one();
    }
#else
    template<typename F, typename... Args>
    void enqueue(F&& f, Args&&... args)
    {
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)();
    }
#endif

    std::vector<Range> splitRange(size_t length,
                                  size_t minSize = 1, size_t jobsPerThread = 1);

    std::vector<Range> splitRange(Range r,
                                  size_t minSize = 1, size_t jobsPerThread = 1);

    void wait();

    int activeThreadCount();

    int queueLength();

    size_t threadCount();

private:
    ThreadPool();

    void threadFunc()
    {
        while (true) {
            std::function<void()> task;

            {
                //std::cout << "trying lock" << std::endl;
                std::unique_lock<std::mutex> lock(m_mutex);
                //std::cout << m_tasks.size() << ' ' << m_workingThreadCount << std::endl;

                if(m_tasks.empty() && m_workingThreadCount == 0)
                {
                    //std::cout << "notified pool done" << std::endl;
                    m_poolDone.notify_all();
                }

                m_cv.wait(lock, [this] {
                    return !m_running || !m_tasks.empty();
                });

                if (!m_running && m_tasks.empty()) {
                    break;
                }

                task = std::move(m_tasks.front());
                m_workingThreadCount++;
                m_tasks.pop();
            }

            //std::cout << "task started" << std::endl;;
            task();
            m_workingThreadCount--;
            //std::cout << "task done" << std::endl;
        }
    }
    std::vector<std::thread> m_threads;
    std::queue<std::function<void()>> m_tasks;

    std::condition_variable m_cv;
    std::condition_variable m_poolDone;
    std::mutex m_mutex;
    std::atomic_bool m_running;
    std::atomic_uint m_workingThreadCount;
};
#endif // THREADPOOL_H
