# -*- coding: utf-8 -*-
import scrapy


class RedditbotSpider(scrapy.Spider):
    name = 'redditbot'
    allowed_domains = ['www.reddit.com/r/gameofthrones/']
    start_urls = ['http://www.reddit.com/r/gameofthrones//']

    def parse(self, response):
        # Extract the content using css selectors
        titles = response.css('.title.may-blank::text').extract()
        votes = response.css('.score.unvoted::text').extract()
        times = response.css('time::attr(title)').extract()
        comments = response.css('.comments::text').extract()

        for item in zip(titles, votes, times, comments):
            scraped_info = {
                'title': item[0],
                'vote': item[1],
                'created at': item[2],
                'comments': item[3],
            }

            yield scraped_info

