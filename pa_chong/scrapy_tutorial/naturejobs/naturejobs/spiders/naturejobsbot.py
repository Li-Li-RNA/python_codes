# -*- coding: utf-8 -*-
import scrapy


class NaturejobsbotSpider(scrapy.Spider):
    name = 'naturejobsbot'
    # allowed_domains = ['https://www.nature.com/naturejobs/science/jobs']
    start_urls = ['https://www.nature.com/naturejobs/science/jobs']

    def parse(self, response):
        # Extract the content using css selectors
        employer = response.css('.employer::text').extract()
        locale = response.css('.locale::text').extract()
        time = response.css('.when::text').extract()

        for item in zip(employer, locale, time):
            scraped_info = {
                'employer': item[0],
                'locale': item[1],
                'created at': item[2]
            }

            yield scraped_info

        next_page = response.css('div.pagination.next a::attr(href)').extract_first()
        if next_page is not None:
            next_page = response.urljoin(next_page)
            yield scrapy.Request(next_page, callback=self.parse)
