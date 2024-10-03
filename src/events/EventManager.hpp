#ifndef EVENT_MANAGER_HPP
#define EVENT_MANAGER_HPP
#include <vector>
#include <unordered_map>

class Event;

class EventManager {
    public:
        EventManager(void);
        void registerEvent(Event* e, int interval);
        void initialize();
        void call(int iteration);
    private:
        std::unordered_map<int, std::vector<Event*>> events;
};

#endif