#ifndef EVENTS_H
#define EVENTS_H

#include <functional>

// Define custom event types or function types as needed
using FileSelectionChangedEvent = std::function<void(const std::string&)>;
using RefreshEvent = std::function<void()>;
using SliderValueChangedEvent = std::function<void(int)>;
using SerializeEvent = std::function<void()>;

#endif // EVENTS_H
