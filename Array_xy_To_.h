#pragma once

//template<typename T>sf::VertexArray Array_xy_To_(std::vector<std::pair<T, T>>& Array_xy) {
//	sf::VertexArray 
//
//	for (auto const& : Array_xy) {
//        for (float x = 0; x < 10.0f; x += 0.1f) {
//            float y = F((float)x);
//            float screenX = x * xScale;
//            float screenY = yOffset - y * amplitude;
//            graph.append(sf::Vertex{ sf::Vector2f(screenX, screenY),sf::Color::Red });
//        }
//	}
//
//}

#if 1
#include <SFML/Graphics.hpp>
#include <cmath>
#include <functional>
template<class T>sf::VertexArray funct(std::function<T(T)>& F) {
    sf::VertexArray graph(sf::PrimitiveType::LineStrip, 0);
    float amplitude = 50.0f;
    float xScale = 80.0f;
    float yOffset = 300.0f;

    for (float x = 0; x < 10.0f; x += 0.1f) {
        float y = F((float)x);
        float screenX = x * xScale;
        float screenY = yOffset - y * amplitude;
        graph.append(sf::Vertex{ sf::Vector2f(screenX, screenY),sf::Color::Red });
    }
    return graph;
}
template<class T>void foo(T f) {
    sf::RenderWindow window(sf::VideoMode({ 800, 600 }), "График функции y = sin(x)");

    while (window.isOpen()) {

        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>() ||
                (event->is<sf::Event::KeyPressed>() &&
                    event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Escape))
            {
                window.close();
            }
        }

        window.clear(sf::Color::White);
        window.draw(funct(static_cast<std::function<double(double)>>(f)));
        /*window.draw(funct(static_cast<double(*)(double)>(std::sin)));*/
        window.display();
    }

}

int main() {

    foo(static_cast<double(*)(double)>(std::cos));

    return 0;
}
#endif