#coding=utf-8
import re
import urllib

def getHtml(url):
    page = urllib.urlopen(url)
    html = page.read()
    return html

html = getHtml("http://www.jianshu.com")
reg = r'<h4 class="title"><a target="_blank" href="(.*?)">(.*?)</a></h4>[\n][\s]+<div class="list-footer">[\n][\s]+<a target="_blank" href="/p/.*?">[\n][\s]+(.*?)[\n]</a>[\s]+<a target="_blank" href=".*?">[\n][\s]+(.*?)[\n]</a>[\s]+<span>(.*?)</span>'
hotre = re.compile(reg)
artlist = re.findall(hotre, html)

for article in artlist:
    for com in article:
        if com.startswith("/p/"):
            print "http://www.jianshu.com"+com
        else:
            print com


